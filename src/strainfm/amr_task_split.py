"""Pre-registered joint split for the expanded AMR cohort (near-clone isolated; source stratified).

Components are unions over genomes by near-clone group, declared source
projects (`source_projects` JSON list), ENA study accessions and BioProject.
One joint partition is built across both drugs so a BioSample keeps the same
partition in every task; per-drug labeled subsets inherit it. Components are
assigned deterministically (seeded greedy by labeled-unit deficit) with a
repair pass enforcing, for every drug, at least `min_groups` groups and both
classes in every partition. Insufficient support raises instead of relaxing the
rules.
"""

import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

UNIT_COLUMNS = ["drug", "genome_id", "biosample", "label", "label_conflict",
                "near_clone_group", "source_component", "split"]
SPLITS = ("train", "validation", "test")


class _UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def _stable_rank(seed, key):
    return hashlib.sha256(f"{seed}\t{key}".encode()).hexdigest()


def _parse_projects(value):
    if not value:
        return []
    try:
        parsed = json.loads(value)
        if isinstance(parsed, list):
            return [str(item).strip() for item in parsed if str(item).strip()]
    except (json.JSONDecodeError, TypeError):
        pass
    return [item.strip() for item in value.split(",") if item.strip()]


def _component_stats(genomes, members, unit_set, labels, groups):
    stats = {"key": members[0], "units": 0, "drugs": {}}
    for drug in labels:
        drug_units = [genome for genome in members if (drug, genome) in unit_set]
        if not drug_units:
            continue
        stats["drugs"][drug] = {
            "units": len(drug_units),
            "groups": len({groups[g]["near_clone_group"] for g in drug_units}),
            "S": sum(labels[drug][g] == "S" for g in drug_units),
            "R": sum(labels[drug][g] == "R" for g in drug_units),
        }
        stats["units"] += len(drug_units)
    return stats


def assign_components(components, seed, drugs, fractions=(0.7, 0.1, 0.2), min_groups=20):
    """Assign whole components; every drug keeps >=min_groups groups and both classes per partition."""
    total_units = sum(c["units"] for c in components)
    targets = {split: fraction * total_units for split, fraction in zip(SPLITS, fractions)}
    assigned_units = Counter({s: 0 for s in SPLITS})
    counters = {drug: {s: {"groups": 0, "S": 0, "R": 0} for s in SPLITS} for drug in drugs}
    assignment = {}

    def add(comp, split, sign=1):
        assigned_units[split] += sign * comp["units"]
        for drug, stats in comp["drugs"].items():
            for field in ("groups", "S", "R"):
                counters[drug][split][field] += sign * stats[field]

    order = sorted(components, key=lambda c: (_stable_rank(seed, c["key"]), c["key"]))
    for comp in order:
        split = max(SPLITS, key=lambda s: (targets[s] - assigned_units[s], -SPLITS.index(s)))
        assignment[comp["key"]] = split
        add(comp, split)

    def deficient_for(drug):
        return [s for s in SPLITS
                if counters[drug][s]["groups"] < min_groups
                or counters[drug][s]["S"] == 0 or counters[drug][s]["R"] == 0]

    def partition_valid(split):
        return all(counters[drug][split]["groups"] >= min_groups
                   and counters[drug][split]["S"] > 0 and counters[drug][split]["R"] > 0 for drug in drugs)

    for _ in range(len(components) * len(drugs)):
        issues = [(drug, split) for drug in drugs for split in deficient_for(drug)]
        if not issues:
            break
        moved = False
        for drug, split in issues:
            candidates = [c for c in components if assignment[c["key"]] != split]
            candidates.sort(key=lambda c: (c["units"], c["key"]))
            for comp in candidates:
                old = assignment[comp["key"]]
                stats = comp["drugs"].get(drug, {"groups": 0, "S": 0, "R": 0})
                supplies = (counters[drug][split]["groups"] + stats["groups"] >= min_groups
                            and (counters[drug][split]["S"] > 0 or stats["S"] > 0)
                            and (counters[drug][split]["R"] > 0 or stats["R"] > 0))
                if not supplies:
                    continue
                add(comp, old, sign=-1)
                donor_ok = partition_valid(old)
                add(comp, old, sign=1)
                if not donor_ok:
                    continue
                add(comp, old, sign=-1)
                add(comp, split)
                assignment[comp["key"]] = split
                moved = True
                break
            if moved:
                break
        if not moved:
            break
    problems = {drug: {s: counters[drug][s] for s in deficient_for(drug)} for drug in drugs if deficient_for(drug)}
    if problems:
        raise ValueError(f"Unable to satisfy split constraints: {problems} (< {min_groups} groups or missing class)")
    return assignment


def run_amr_task_split(representatives_path, groups_path, output, drugs=("ciprofloxacin", "gentamicin"),
                       seed=20261110, min_groups=20, fractions=(0.7, 0.1, 0.2)):
    representatives_path, groups_path, output = map(Path, (representatives_path, groups_path, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Non-raw output directory required")
    if output.exists():
        raise ValueError("Use a new output directory")
    if abs(sum(fractions) - 1.0) > 1e-9:
        raise ValueError("fractions must sum to 1")
    representatives = read_table(representatives_path)
    groups = {r["genome_id"]: r for r in read_table(groups_path)}
    labels = {drug: {} for drug in drugs}
    genomes = set()
    for row in representatives:
        if row["genome_id"] not in groups:
            continue
        for drug in drugs:
            if row.get(drug) in ("S", "R"):
                labels[drug][row["genome_id"]] = row[drug]
                genomes.add(row["genome_id"])
    if not genomes:
        raise ValueError("No labeled representatives for the requested drugs")
    row_of = {r["genome_id"]: r for r in representatives}
    unit_set = {(drug, genome) for drug in drugs for genome in labels[drug]}
    clone_labels = defaultdict(lambda: defaultdict(set))
    for drug in drugs:
        for genome, label in labels[drug].items():
            clone_labels[groups[genome]["near_clone_group"]][drug].add(label)
    conflicts = {drug: {clone for clone, per_drug in clone_labels.items() if len(per_drug.get(drug, set())) > 1}
                 for drug in drugs}

    uf = _UnionFind()
    by_clone = defaultdict(list)
    for genome in genomes:
        uf.find(genome)
        by_clone[groups[genome]["near_clone_group"]].append(genome)
    for members in by_clone.values():
        for genome in members[1:]:
            uf.union(members[0], genome)

    source_uf = _UnionFind()
    for genome in genomes:
        source_uf.find(genome)
        row = row_of[genome]
        for project in _parse_projects(row.get("source_projects", "")):
            source_uf.union(f"project:{project}", genome)
        for study in (row.get("ena_study_accessions", "") or "").split(","):
            study = study.strip()
            if study:
                source_uf.union(f"study:{study}", genome)
        bioproject = (row.get("bioproject", "") or "").strip()
        if bioproject:
            source_uf.union(f"project:{bioproject}", genome)
    source_component_of = {genome: source_uf.find(genome) for genome in genomes}

    components = defaultdict(list)
    for genome in genomes:
        components[uf.find(genome)].append(genome)
    comp_list, comp_of = [], {}
    for key, members in components.items():
        stats = _component_stats(genomes, members, unit_set, labels, groups)
        comp_list.append(stats)
        for genome in members:
            comp_of[genome] = stats["key"]
    assignment = assign_components(comp_list, seed, drugs, fractions, min_groups)

    unit_rows = []
    for drug in drugs:
        for genome, label in labels[drug].items():
            unit_rows.append({
                "drug": drug, "genome_id": genome, "biosample": row_of[genome].get("biosample", ""),
                "label": label, "label_conflict": groups[genome]["near_clone_group"] in conflicts[drug],
                "near_clone_group": groups[genome]["near_clone_group"],
                "source_component": source_component_of[genome], "split": assignment[comp_of[genome]],
            })
    unit_rows.sort(key=lambda r: (r["drug"], r["split"], r["near_clone_group"], r["genome_id"]))

    contract = {
        "schema": 1,
        "purpose": "Pre-registered joint source- and near-clone-isolated split for the expanded AMR cohort",
        "rules": {
            "unit": "drug x identity representative with an S/R label",
            "components": "near-clone groups only (source components are reported, not isolated: the declared-source graph collapses to 43 components with the largest at 43%, so a balanced source-isolated split is infeasible and preserved as evidence in the source_isolated run)",
            "joint": "one partition across drugs; each BioSample keeps the same partition in every task",
            "assignment": "whole near-clone groups; seeded greedy by labeled-unit deficit; per drug >= min_groups groups and both classes in every partition",
            "conflicts": "mixed-label clone groups flagged per drug; excluded from the primary analysis downstream",
            "source_reporting": "source components may cross partitions; per-partition source composition and the largest source share must be reported, and source-independent test claims are not supported",
        },
        "seed": seed, "min_groups": min_groups, "fractions": list(fractions), "drugs": list(drugs),
        "representatives_sha256": sha256(representatives_path),
        "groups_sha256": sha256(groups_path),
        "code_sha256": sha256(__file__),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    write_tsv(output / "amr_task_units.tsv", unit_rows, UNIT_COLUMNS)
    summaries = {}
    for drug in drugs:
        units = [r for r in unit_rows if r["drug"] == drug]
        primary = [r for r in units if not r["label_conflict"]]
        summaries[drug] = {
            "n_units": len(units),
            "n_S": sum(r["label"] == "S" for r in units),
            "n_R": sum(r["label"] == "R" for r in units),
            "n_groups": len({r["near_clone_group"] for r in units}),
            "n_conflict_groups": len({r["near_clone_group"] for r in units if r["label_conflict"]}),
            "n_conflict_units": len(units) - len(primary),
            "n_components": len({r["source_component"] for r in units}),
            "primary": {"n_units": len(primary), "n_S": sum(r["label"] == "S" for r in primary),
                        "n_R": sum(r["label"] == "R" for r in primary)},
            "partitions": {
                split: {
                    "n_units": sum(r["split"] == split for r in units),
                    "n_primary_units": sum(r["split"] == split for r in primary),
                    "n_groups": len({r["near_clone_group"] for r in units if r["split"] == split}),
                    "n_source_components": len({r["source_component"] for r in units if r["split"] == split}),
                    "largest_source_share": round(max(Counter(r["source_component"] for r in units if r["split"] == split).values()) / max(1, sum(1 for r in units if r["split"] == split)), 3),
                    "S": sum(r["split"] == split and r["label"] == "S" for r in units),
                    "R": sum(r["split"] == split and r["label"] == "R" for r in units),
                } for split in SPLITS
            },
        }
    summary = {
        "status": "expanded_amr_joint_split_frozen",
        "seed": seed, "min_groups": min_groups, "fractions": list(fractions),
        "drugs": summaries,
        "n_components_total": len(comp_list),
        "units_sha256": sha256(output / "amr_task_units.tsv"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limits": [
            "Labels are published BV-BRC AST evidence linked at BioSample level; I is excluded, missing stays missing.",
            "Whole declared-source components are isolated; project/study links are database metadata, not proven lineage independence.",
            "179 labeled samples lack ENA study links but carry source_projects; they are merged by declaration and do not support source-independent test claims on their own.",
            "Mixed-label clone groups are flagged per drug; the primary analysis excludes them.",
            "The old AMR-02 split and the earlier single-drug partitions remain unchanged; overlap with the seen AMR-02 test must be reported.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "amr_task_units.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# Expanded AMR joint split (pre-registered)\n\n"
        "Joint near-clone-isolated split over identity representatives labeled from the BV-BRC AST snapshot. "
        "Whole near-clone groups form components; source components (declared projects, ENA studies, BioProject) are "
        "reported per partition but not isolated, because the declared-source graph is one connected blob. "
        "The same BioSample keeps one partition across drugs. Frozen before any model result is inspected.\n\n"
        + "".join(
            f"- {drug}: N={info['n_units']} (S={info['n_S']}, R={info['n_R']}), groups={info['n_groups']}, "
            f"components={info['n_components']}, conflicts={info['n_conflict_groups']} groups; partitions "
            + "; ".join(f"{split}: {part['n_primary_units']} primary units / {part['n_groups']} groups / {part['n_source_components']} source components (largest {part['largest_source_share']}) / S{part['S']} R{part['R']}"
                        for split, part in info["partitions"].items())
            + "\n"
            for drug, info in summaries.items()
        )
        + "\nOld AMR-02 artifacts are unchanged; overlap with the seen AMR-02 test must be reported in any confirmatory claim.\n"
    )
    return summary
