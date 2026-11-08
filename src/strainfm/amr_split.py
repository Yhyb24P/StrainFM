"""Pre-registered AMR task split over labeled E. coli BioSamples (D0-27/D0-28).

Rules frozen before any model result is inspected:
- one unit per (drug, BioSample, versioned assembly); paired-accession labels are
  the primary evidence, direct file-accession matching is recorded as the
  sensitivity subset;
- the grouping unit is the near-clone group; whole groups stay in one partition;
- source components may cross partitions (documented limitation);
- clone groups with mixed S/R labels are flagged and excluded from the primary
  analysis (sensitivity may include them);
- partitions are 70/10/20 by group count with a fixed seed and at least
  `min_groups` groups per partition; insufficient support raises instead of
  silently relaxing the rule.
"""

import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .snapshot import atomic_json, prepare_output

UNIT_COLUMNS = [
    "drug", "genome_id", "biosample", "label", "label_source",
    "near_clone_group", "source_group", "label_conflict", "split",
]
DIRECT_LABEL_SOURCES = ("file_accession",)


def _stable_rank(seed, drug, group):
    return hashlib.sha256(f"{seed}\t{drug}\t{group}".encode()).hexdigest()


def assign_groups(group_labels, seed, drug, fractions=(0.7, 0.1, 0.2), min_groups=20):
    """Deterministically assign whole groups to train/validation/test.

    group_labels: {group: sorted list of labels in that group} (labels are "S"/"R").
    Greedy by seeded order: each group goes to the partition with the largest
    remaining deficit (by group count). Returns {group: split}.
    """
    splits = ("train", "validation", "test")
    total = len(group_labels)
    targets = {split: fraction * total for split, fraction in zip(splits, fractions)}
    counts = Counter({split: 0 for split in splits})
    assignment = {}
    ordered = sorted(group_labels, key=lambda group: (_stable_rank(seed, drug, group), group))
    for group in ordered:
        split = max(splits, key=lambda name: (targets[name] - counts[name], -splits.index(name)))
        assignment[group] = split
        counts[split] += 1
    short = {split: counts[split] for split in splits if counts[split] < min_groups}
    if short:
        raise ValueError(f"Insufficient group support for {drug}: {short} (< {min_groups})")
    return assignment


def run_amr_split(labels_path, cohort_manifest_path, split_manifest_path, nonrepresentatives_path,
                  output, drugs=("gentamicin", "ciprofloxacin"), seed=20261108, min_groups=20):
    labels_path, cohort_manifest_path, split_manifest_path, nonrepresentatives_path, output = map(
        Path, (labels_path, cohort_manifest_path, split_manifest_path, nonrepresentatives_path, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Non-raw output directory required")
    if output.exists():
        raise ValueError("Use a new output directory")
    cohort = read_table(cohort_manifest_path)
    split_rows = read_table(split_manifest_path)
    nonrepresentatives = read_table(nonrepresentatives_path)
    clone_of = {r["genome_id"]: r["near_clone_group"] for r in split_rows}
    clone_of.update({r["genome_id"]: r["near_clone_group"] for r in nonrepresentatives})
    source_of_clone = {r["near_clone_group"]: r.get("source_group", "") for r in split_rows}
    accession_owner = {}
    for row in cohort:
        for field in ("file_accession", "current_accession", "paired_accession"):
            value = row.get(field, "")
            if value and value not in ("missing", "unknown") and value not in accession_owner:
                accession_owner[value] = row
    labels = {}
    with labels_path.open(newline="") as handle:
        import csv
        reader = csv.reader(handle)
        header = next(reader)
        drugs = list(drugs)
        indices = {drug: header.index(drug) - 1 for drug in drugs}
        for row in reader:
            labels[row[0]] = row[1:]
    matched = [accession_owner[acc] for acc in labels if acc in accession_owner]
    eligible = [row for row in matched if row.get("eligible") != "False"]
    unit_rows, summaries = [], {}
    for drug in drugs:
        index = indices[drug]
        units = []
        for row in eligible:
            label, source = None, None
            for field in ("file_accession", "current_accession", "paired_accession"):
                accession = row.get(field, "")
                if accession in labels:
                    value = labels[accession][index]
                    if value in ("0", "0.0", "1", "1.0"):
                        label = "R" if value in ("1", "1.0") else "S"
                        source = field
                        break
            if label is None:
                continue
            genome_id = row["file_accession"]
            group = clone_of.get(genome_id) or clone_of.get(row.get("genome_id", ""))
            if group is None:
                raise ValueError(f"Labeled assembly without near-clone group: {genome_id}")
            units.append({"genome_id": genome_id, "biosample": row.get("biosample", ""),
                          "label": label, "label_source": source,
                          "near_clone_group": group,
                          "source_group": source_of_clone.get(group, "")})
        groups = defaultdict(list)
        for unit in units:
            groups[unit["near_clone_group"]].append(unit["label"])
        conflicts = {group for group, values in groups.items() if len(set(values)) > 1}
        assignment = assign_groups(groups, seed, drug, min_groups=min_groups)
        for unit in units:
            unit["label_conflict"] = unit["near_clone_group"] in conflicts
            unit["split"] = assignment[unit["near_clone_group"]]
            unit["drug"] = drug
        unit_rows.extend(units)
        primary = [u for u in units if not u["label_conflict"]]
        summaries[drug] = {
            "n_units": len(units),
            "n_S": sum(u["label"] == "S" for u in units),
            "n_R": sum(u["label"] == "R" for u in units),
            "n_direct_file_accession": sum(u["label_source"] in DIRECT_LABEL_SOURCES for u in units),
            "n_groups": len(groups),
            "n_conflict_groups": len(conflicts),
            "n_conflict_units": len(units) - len(primary),
            "primary_n_units": len(primary),
            "primary_n_S": sum(u["label"] == "S" for u in primary),
            "primary_n_R": sum(u["label"] == "R" for u in primary),
            "partitions": {
                split: {
                    "n_units": sum(u["split"] == split for u in units),
                    "n_primary_units": sum(u["split"] == split for u in units if not u["label_conflict"]),
                    "n_groups": sum(1 for group, name in assignment.items() if name == split),
                    "n_source_components": len({source_of_clone.get(group, "") for group, name in assignment.items() if name == split}),
                } for split in ("train", "validation", "test")
            },
        }
    unit_rows.sort(key=lambda row: (row["drug"], row["near_clone_group"], row["genome_id"]))
    contract = {
        "schema": 1,
        "purpose": "Pre-registered near-clone-group split for AMR task evaluation; frozen before model results",
        "rules": {
            "unit": "drug x BioSample x versioned assembly; one row per unit",
            "group": "near_clone_group; whole group in one partition",
            "source_isolation": False,
            "source_note": "source components may cross partitions by user decision (D0-28); documented in limitations",
            "conflicts": "mixed S/R within a clone group flagged; primary analysis excludes them",
            "partitioning": "70/10/20 by group count, deterministic seeded greedy; min_groups enforced",
        },
        "seed": seed,
        "min_groups": min_groups,
        "drugs": list(drugs),
        "labels_sha256": sha256(labels_path),
        "cohort_manifest_sha256": sha256(cohort_manifest_path),
        "split_manifest_sha256": sha256(split_manifest_path),
        "nonrepresentatives_sha256": sha256(nonrepresentatives_path),
        "code_sha256": sha256(__file__),
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    write_tsv(output / "amr_split_units.tsv", unit_rows, UNIT_COLUMNS)
    summary = {
        "status": "amr_split_frozen",
        "seed": seed,
        "min_groups": min_groups,
        "drugs": summaries,
        "n_quality_excluded_matched": len(matched) - len(eligible),
        "n_matched_assemblies": len(matched),
        "units_sha256": sha256(output / "amr_split_units.tsv"),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limits": [
            "Labels are published binary AST results (paired-accession evidence); direct-only is a sensitivity subset (496 per drug).",
            "Source components cross partitions by user decision; shared project is not lineage proof.",
            "Mixed S/R clone groups are excluded from the primary analysis and reported.",
            "The split is frozen before model results; the global retrieval split is unchanged.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "amr_split_units.tsv", "contract.json")})
    (output / "README.md").write_text(
        "# AMR task split (pre-registered)\n\n"
        "Near-clone-group-isolated split for the two drugs that pass the 500/50 count gate under "
        "paired-accession evidence. Frozen before any frozen-representation result is inspected.\n\n"
        + "".join(
            f"- {drug}: N={info['n_units']} (S={info['n_S']}, R={info['n_R']}), groups={info['n_groups']}, "
            f"conflict groups={info['n_conflict_groups']}; partitions "
            + "; ".join(f"{split}: {part['n_primary_units']} primary units / {part['n_groups']} groups"
                        for split, part in info["partitions"].items())
            + "\n"
            for drug, info in summaries.items()
        )
        + "\nThe global retrieval split is unchanged. Source components may cross partitions (D0-28).\n"
    )
    return summary
