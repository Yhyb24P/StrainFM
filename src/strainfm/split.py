"""Freeze source-isolated splits after exact identity and near-clone grouping."""

import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import UnionFind, rank_quality, read_table
from .metadata import known
from .snapshot import atomic_json


def run_split(representatives_path, groups_path, cohort_manifest_path, output, seed=20261105):
    representatives_path, groups_path, cohort_manifest_path, output = map(Path, (representatives_path, groups_path, cohort_manifest_path, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    rows = read_table(representatives_path)
    groups = {r["genome_id"]: r for r in read_table(groups_path)}
    if set(groups) != {r["genome_id"] for r in rows}:
        raise ValueError("Grouping inventory differs from cohort representatives")
    clones = defaultdict(list)
    for row in rows:
        row.update(groups[row["genome_id"]])
        clones[row["near_clone_group"]].append(row)
    identity_members = defaultdict(list)
    for member in read_table(cohort_manifest_path):
        if member["eligible"] == "True":
            identity_members[member["identity_representative"]].append(member)
    if set(identity_members) != set(groups):
        raise ValueError("Full identity membership differs from grouping inventory")
    retained = []
    project_union = UnionFind({member["bioproject"] for members in identity_members.values() for member in members if known(member.get("bioproject", ""))})
    for component in clones.values():
        representative = min(component, key=rank_quality)
        members = [member for row in component for member in identity_members[row["genome_id"]]]
        sources = sorted({member["bioproject"] for member in members if known(member.get("bioproject", ""))})
        representative["source_projects"] = json.dumps(sources, separators=(",", ":"))
        representative["identity_and_clone_member_count"] = len(members)
        for source in sources[1:]:
            project_union.union(sources[0], source)
        retained.append(representative)
    projects, unknown = defaultdict(list), []
    for row in retained:
        sources = json.loads(row["source_projects"])
        if sources:
            projects[project_union.find(sources[0])].append(row)
        else:
            unknown.append(row)
    if len(projects) < 3:
        raise ValueError("Fewer than three source components after preserving all identity/clone member projects; independent splits unsupported")
    fractions = {"train": 0.8, "validation": 0.1, "test": 0.1}
    counts = Counter({"train": len(unknown), "validation": 0, "test": 0})
    total = len(retained)
    targets = {split: fraction * total for split, fraction in fractions.items()}
    def stable_rank(project):
        return hashlib.sha256(f"{seed}\t{project}".encode()).hexdigest()
    for project in sorted(projects, key=lambda p: (-len(projects[p]), stable_rank(p))):
        split = max(fractions, key=lambda s: targets[s] - counts[s])
        for row in projects[project]:
            row["split"] = split
            row["source_group"] = project
        counts[split] += len(projects[project])
    for row in unknown:
        row.update(split="train", source_group="unknown")
    for row in retained:
        row.update(pretraining_overlap="unknown", overlap_note="No complete upstream accession membership list verified")
    # Verify actual emitted groups, rather than relying on split algorithm names.
    for field in ("genome_id", "file_accession", "biosample", "identity_group", "near_clone_group", "bioproject"):
        owners = {}
        for row in retained:
            value = row.get(field, "")
            if known(value):
                if value in owners and owners[value] != row["split"]:
                    raise ValueError(f"Cross-split leakage in {field}: {value}")
                owners[value] = row["split"]
    source_owners = {}
    for row in retained:
        for project in json.loads(row["source_projects"]):
            if project in source_owners and source_owners[project] != row["split"]:
                raise ValueError(f"Cross-split source leakage in merged members: {project}")
            source_owners[project] = row["split"]
    if not all(counts[s] for s in fractions):
        raise ValueError("Empty evaluation split")
    training_neighborhoods = Counter(r["sequence_neighborhood"] for r in retained if r["split"] == "train")
    queries = [{"genome_id": row["genome_id"], "split": row["split"], "positive_definition": "same_Mash_sequence_neighborhood", "gallery_split": "train", "gallery_positive_count": training_neighborhoods[row["sequence_neighborhood"]], "eligible_for_positive_supported_retrieval": training_neighborhoods[row["sequence_neighborhood"]] > 0} for row in retained if row["split"] != "train"]
    supported = Counter(r["split"] for r in queries if r["eligible_for_positive_supported_retrieval"])
    if not supported["validation"] or not supported["test"]:
        raise ValueError("No positive-supported retrieval queries in validation or test")
    output.mkdir(parents=True)
    columns = list(dict.fromkeys(key for row in retained for key in row))
    write_tsv(output / "split_manifest.tsv", sorted(retained, key=lambda r: r["genome_id"]), columns)
    write_tsv(output / "retrieval_query_support.tsv", queries, list(queries[0]))
    kept = {r["genome_id"] for r in retained}
    removed = [{"genome_id": row["genome_id"], "near_clone_group": row["near_clone_group"], "reason": "nonrepresentative_near_clone"} for row in rows if row["genome_id"] not in kept]
    write_tsv(output / "near_clone_nonrepresentatives.tsv", removed, ["genome_id", "near_clone_group", "reason"])
    source_counts = {split: sum(owner == split for owner in source_owners.values()) for split in fractions}
    neighborhoods = {split: {r["sequence_neighborhood"] for r in retained if r["split"] == split} for split in fractions}
    summary = {"status": "candidate_splits_frozen_pending_target_species_confirmation", "seed": seed, "identity_representatives": len(rows), "near_clone_representatives": len(retained), "split_counts": dict(counts), "source_project_counts": source_counts, "source_unknown_train_only": len(unknown), "retrieval_positive_supported_queries": dict(supported), "retrieval_queries_without_gallery_positive": {s: counts[s] - supported[s] for s in ("validation", "test")}, "sequence_neighborhood_counts": {s: len(v) for s, v in neighborhoods.items()}, "test_neighborhoods_seen_in_train": len(neighborhoods["test"] & neighborhoods["train"]), "test_neighborhoods_absent_from_train": len(neighborhoods["test"] - neighborhoods["train"]), "verified_no_cross_split_identity_nearclone_or_known_project": True, "input_sha256": {str(p): sha256(p) for p in (representatives_path, groups_path)}, "code_sha256": sha256(__file__), "split_manifest_sha256": sha256(output / "split_manifest.tsv"), "limitations": ["Source-isolated within-species retrieval is not automatically unseen-lineage evaluation.", "Neighborhood positives derive from Mash grouping and diagnose genetic structure, not independent phenotypes.", "All upstream pretraining exposure remains unknown; report this evaluation stratum explicitly.", "Approximate MinHash near-clone groups are operational, not proof of epidemiological independence.", "Fractions are targets; whole-project allocation prioritizes source isolation."]}
    summary["input_sha256"][str(cohort_manifest_path)] = sha256(cohort_manifest_path)
    summary.update(source_components=len(projects), largest_source_component_representatives=max(len(v) for v in projects.values()), full_member_source_constraints_verified=True)
    atomic_json(output / "summary.json", summary)
    return summary
