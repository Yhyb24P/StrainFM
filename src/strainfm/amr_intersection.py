"""Intersect the BacBench AMR binary labels with the frozen E. coli cohort before any modelling."""

import csv
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .cohort import read_table
from .metadata import known
from .snapshot import atomic_json, prepare_output


def _read_labels(path):
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        genomes, drugs = [], header[1:]
        for row in reader:
            if not row:
                continue
            values = row[1:]
            normalized = {}
            for drug, value in zip(drugs, values):
                value = value.strip()
                normalized[drug] = "0" if value in ("0", "0.0") else "1" if value in ("1", "1.0") else value
            genomes.append((row[0].strip(), normalized))
    return genomes, drugs


def _base_accession(accession):
    return accession.rsplit(".", 1)[0] if accession.count(".") == 1 else accession


def run_amr_intersection(labels_path, cohort_manifest_path, split_manifest_path, output):
    labels_path, cohort_manifest_path, split_manifest_path, output = map(Path, (labels_path, cohort_manifest_path, split_manifest_path, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    label_rows, drugs = _read_labels(labels_path)
    cohort = read_table(cohort_manifest_path)
    representatives = {r["genome_id"]: r for r in read_table(split_manifest_path)}
    accession_owner = {}
    for row in cohort:
        for field in ("file_accession", "current_accession", "paired_accession"):
            value = row.get(field, "")
            if known(value) and value not in accession_owner:
                accession_owner[value] = (row["file_accession"], field)
    base_owner = defaultdict(list)
    for accession in accession_owner:
        base_owner[_base_accession(accession)].append(accession)
    exact_hits, paired_hits, version_mismatch = {}, {}, []
    unmatched = 0
    for accession, values in label_rows:
        if accession in accession_owner:
            owner, field = accession_owner[accession]
            exact_hits[owner] = ("direct" if field == "file_accession" else f"via_{field}", values)
        elif _base_accession(accession) in base_owner:
            candidates = base_owner[_base_accession(accession)]
            version_mismatch.append({"label_accession": accession, "cohort_accessions": ";".join(sorted(candidates))})
        else:
            unmatched += 1
    matched_cohort_rows = {row["file_accession"]: row for row in cohort if row["file_accession"] in exact_hits}
    for row in matched_cohort_rows.values():
        if "paired_accession" in exact_hits[row["file_accession"]][0]:
            paired_hits[row["file_accession"]] = exact_hits[row["file_accession"]]
    matched_representatives = {gid: representatives[gid] for gid in representatives if gid in matched_cohort_rows}
    biosample_groups = defaultdict(set)
    for gid, row in matched_representatives.items():
        biosample = row.get("biosample", "")
        if known(biosample):
            biosample_groups[biosample].add(gid)
    duplicate_biosamples = {bs: sorted(gids) for bs, gids in biosample_groups.items() if len(gids) > 1}
    matched_rows = []
    drug_stats = {drug: Counter() for drug in drugs}
    for gid, row in sorted(matched_representatives.items()):
        match_type, values = exact_hits[gid]
        labeled = {drug: value for drug, value in values.items() if value in ("0", "1")}
        split = row["split"]
        for drug, value in labeled.items():
            drug_stats[drug][f"n_{ 'R' if value == '1' else 'S'}"] += 1
            drug_stats[drug][f"{split}_{'R' if value == '1' else 'S'}"] += 1
            drug_stats[drug][f"{split}_total"] += 1
            drug_stats[drug]["n_labeled"] += 1
        matched_rows.append({
            "genome_id": gid, "split": split, "biosample": row.get("biosample", ""),
            "bioproject": row.get("bioproject", ""), "assembly_level": row.get("assembly_level", ""),
            "near_clone_group": row.get("near_clone_group", ""), "sequence_neighborhood": row.get("sequence_neighborhood", ""),
            "match_type": match_type, "n_drugs_labeled": len(labeled),
        })
    matched_rows.sort(key=lambda r: r["genome_id"])
    drug_rows = []
    for drug in drugs:
        stats = drug_stats[drug]
        n_s, n_r = stats["n_S"], stats["n_R"]
        drug_rows.append({
            "drug": drug, "n_labeled": stats["n_labeled"], "n_S": n_s, "n_R": n_r,
            "train_S": stats["train_S"], "train_R": stats["train_R"], "train_total": stats["train_total"],
            "validation_S": stats["validation_S"], "validation_R": stats["validation_R"], "validation_total": stats["validation_total"],
            "test_S": stats["test_S"], "test_R": stats["test_R"], "test_total": stats["test_total"],
            "paper_rule_min_500_and_50_per_class": stats["n_labeled"] >= 500 and min(n_s, n_r) >= 50,
        })
    drug_rows.sort(key=lambda r: -r["n_labeled"])
    funnel = [
        {"stage": "cohort_records", "n": len(cohort)},
        {"stage": "label_rows", "n": len(label_rows)},
        {"stage": "exact_accession_matches_cohort", "n": len(matched_cohort_rows)},
        {"stage": "version_mismatch_only", "n": len(version_mismatch)},
        {"stage": "unmatched_label_rows", "n": unmatched},
        {"stage": "matched_split_representatives", "n": len(matched_representatives)},
        {"stage": "unique_biosamples_matched", "n": len(biosample_groups)},
        {"stage": "biosamples_with_multiple_representatives", "n": len(duplicate_biosamples)},
    ]
    contract = {
        "schema": 1,
        "purpose": "Pre-modelling intersection of BacBench AMR binary labels with the frozen E. coli cohort; no model inference",
        "labels_sha256": sha256(labels_path),
        "cohort_manifest_sha256": sha256(cohort_manifest_path),
        "split_manifest_sha256": sha256(split_manifest_path),
        "code_sha256": sha256(__file__),
        "label_direction": "0=susceptible, 1=resistant as published (AST S/R only, ambiguous excluded); direction not re-verified against source AST in this audit",
        "matching_rule": "Exact versioned accession against file/current/paired accession; version-stripped matches counted separately as diagnostics and never merged",
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    write_tsv(output / "funnel.tsv", funnel, ["stage", "n"])
    write_tsv(output / "matched_representatives.tsv", matched_rows, list(matched_rows[0]) if matched_rows else ["genome_id"])
    write_tsv(output / "drug_counts_representatives.tsv", drug_rows, list(drug_rows[0]))
    write_tsv(output / "version_mismatch_matches.tsv", version_mismatch, ["label_accession", "cohort_accessions"])
    summary = {
        "status": "amr_intersection_complete",
        "funnel": {r["stage"]: r["n"] for r in funnel},
        "matched_representatives_by_split": dict(Counter(r["split"] for r in matched_rows)),
        "n_drugs": len(drugs),
        "n_drugs_with_any_label": sum(r["n_labeled"] > 0 for r in drug_rows),
        "n_drugs_passing_paper_rule": sum(r["paper_rule_min_500_and_50_per_class"] for r in drug_rows),
        "n_duplicate_biosample_groups": len(duplicate_biosamples),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limitations": [
            "Intersection counts are label availability, not verified AST quality or independent-sample counts.",
            "GCA/GCF matches rely on NCBI paired accession metadata, not on sequence identity between the two records.",
            "Same BioSample with multiple assemblies is reported, not silently collapsed.",
            "Binary labels exclude ambiguous (I) results by the published dataset construction; per-drug conflicts are not observable in this file.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "funnel.tsv", "matched_representatives.tsv", "drug_counts_representatives.tsv", "version_mismatch_matches.tsv", "contract.json")})
    return summary
