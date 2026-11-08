"""Trace AMR label loss and count task feasibility; never select or split a training set."""

import csv
import json
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .audit import safe_path, sha256, write_tsv
from .cohort import read_table
from .metadata import known
from .snapshot import atomic_json, prepare_output

ACCESSION_FIELDS = ("file_accession", "current_accession", "paired_accession")
SPLITS = ("train", "validation", "test")


def _read_labels(path):
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if header[0] != "genome_name" or len(set(header)) != len(header) or len(header) < 2:
            raise ValueError("Expected unique drug columns after genome_name")
        genomes, drugs = [], header[1:]
        for number, row in enumerate(reader, 2):
            if not row:
                continue
            if len(row) != len(header) or not row[0].strip():
                raise ValueError(f"Malformed label row {number}")
            normalized = {}
            for drug, value in zip(drugs, row[1:]):
                value = value.strip()
                normalized[drug] = "0" if value in ("0", "0.0") else "1" if value in ("1", "1.0") else "" if value.lower() in ("", "nan", "na", "none") else value
            genomes.append((row[0].strip(), normalized))
    return genomes, drugs


def _base_accession(accession):
    return accession.rsplit(".", 1)[0] if accession.count(".") == 1 else accession


def _unique(rows, key):
    indexed = {r[key]: r for r in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"Duplicate {key} in inventory")
    return indexed


def _file_evidence(item):
    row, root = item
    path = safe_path(root, row["target_relpath"])
    try:
        before = path.stat()
        digest = sha256(path)
        after = path.stat()
        stable = (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
        status = "snapshot_content_verified" if stable and digest == row["file_sha256"] else "changed_since_snapshot"
        return status, digest
    except OSError as exc:
        return ("missing" if isinstance(exc, FileNotFoundError) else "unreadable"), ""


def run_amr_intersection(labels_path, cohort_manifest_path, split_manifest_path, output, *, groups_path, genome_root):
    labels_path, cohort_manifest_path, split_manifest_path, output, groups_path, genome_root = map(
        Path, (labels_path, cohort_manifest_path, split_manifest_path, output, groups_path, genome_root))
    if output.exists() or "raw" in output.resolve().parts or output.resolve().is_relative_to(genome_root.resolve()):
        raise ValueError("Use a new output directory outside raw/genome directories")
    inputs = (labels_path, cohort_manifest_path, split_manifest_path, groups_path)
    input_hashes = {str(p): sha256(p) for p in inputs}
    label_rows, drugs = _read_labels(labels_path)
    cohort = _unique(read_table(cohort_manifest_path), "genome_id")
    representatives = _unique(read_table(split_manifest_path), "genome_id")
    groups = _unique(read_table(groups_path), "genome_id")
    clone_representatives = _unique(list(representatives.values()), "near_clone_group")
    if any(r["split"] not in SPLITS for r in representatives.values()):
        raise ValueError("Unknown existing split")
    by_accession = defaultdict(list)
    for accession, values in label_rows:
        by_accession[accession].append(values)
    cohort_accessions = {r.get(f, "") for r in cohort.values() for f in ACCESSION_FIELDS if known(r.get(f, ""))}
    cohort_bases = {_base_accession(a) for a in cohort_accessions}
    version_only = sorted(a for a in by_accession if a not in cohort_accessions and _base_accession(a) in cohort_bases)
    matched = []
    calls = {}
    for gid, row in sorted(cohort.items()):
        hits = {row[f] for f in ACCESSION_FIELDS if row.get(f) in by_accession}
        if not hits:
            continue
        calls[gid] = {d: {v[d] for a in hits for v in by_accession[a] if v[d] != ""} for d in drugs}
        if row["eligible"] == "True":
            identity_rep = row["identity_representative"]
            if identity_rep not in groups:
                raise ValueError(f"Missing sequence group for eligible assembly {gid}")
            near = groups[identity_rep]["near_clone_group"]
            if near not in clone_representatives:
                raise ValueError(f"Missing near-clone representative for {gid}")
            rep = clone_representatives[near]
            if gid in representatives:
                disposition = "retained_split_representative"
            elif row["retained_identity_representative"] != "True":
                disposition = "identity_nonrepresentative"
            else:
                disposition = "near_clone_nonrepresentative"
        else:
            near, rep, disposition = "", {}, "cohort_quality_excluded"
        matched.append({
            "genome_id": gid, "file_accession": row["file_accession"], "biosample": row.get("biosample", ""),
            "bioproject": row.get("bioproject", ""), "label_accessions": json.dumps(sorted(hits)),
            "match_fields": ";".join(f for f in ACCESSION_FIELDS if row.get(f) in hits),
            "label_row_count": sum(len(by_accession[a]) for a in hits),
            "disposition": disposition, "exclusion_reason": row["exclusion_reason"],
            "snapshot_status": row["status"], "eligible": row["eligible"],
            "identity_group": row.get("identity_group", ""), "identity_representative": row.get("identity_representative", ""),
            "near_clone_group": near, "selected_clone_representative": rep.get("genome_id", ""),
            "selected_representative_has_label": rep.get("genome_id") in calls or any(rep.get(f) in by_accession for f in ACCESSION_FIELDS),
            "inherited_split": rep.get("split", ""), "source_group": rep.get("source_group", ""),
            "n_binary_drugs": sum(states in ({"0"}, {"1"}) for states in calls[gid].values()),
            "n_conflicting_drugs": sum({"0", "1"} <= states for states in calls[gid].values()),
            "file_sha256": row["file_sha256"], "target_relpath": row["target_relpath"],
        })
    if not matched:
        raise ValueError("No exact versioned accession matches")
    matched_file_accessions = {r["file_accession"] for r in matched}
    for row in matched:
        row["loss_mechanism"] = row["disposition"]
        if row["disposition"] == "retained_split_representative":
            row["loss_mechanism"] = "legacy_recognized_representative" if row["genome_id"] in matched_file_accessions else "legacy_identifier_key_mismatch"
    contract = {
        "schema": 2, "purpose": "AMR-01 label-loss and per-drug feasibility audit; no training, selection or new split",
        "inputs_sha256": input_hashes, "genome_root": str(genome_root),
        "code_sha256": {p.name: sha256(p) for p in (Path(__file__), Path(__file__).with_name("audit.py"), Path(__file__).with_name("cohort.py"), Path(__file__).with_name("metadata.py"))},
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "label_direction": "Published binary CSV: 0=S, 1=R; underlying AST quality/conditions not reverified",
        "matching_rule": "Exact versioned file/current/paired accession; all matches retained, no label transfer to clone representatives",
        "quality_rule": "Reuse frozen cohort eligibility only when current compressed SHA-256 matches its fully decoded snapshot",
        "count_rule": "Known BioSample, nonconflicting binary label directly on at least one quality-eligible assembly; no replacement assembly selected",
        "paper_count_gate": {"total": 500, "per_class": 50},
        "group_support_rule": "At least 3 groups carrying each class is necessary, not sufficient, for 3-way class coverage; no new split constructed",
        "started_utc": datetime.now(timezone.utc).isoformat(),
    }
    prepare_output(output, contract)
    with ThreadPoolExecutor(max_workers=8) as pool:
        evidence = list(pool.map(_file_evidence, [(cohort[r["genome_id"]], genome_root) for r in matched]))
    for row, (status, digest) in zip(matched, evidence):
        row.update(current_file_status=status, observed_file_sha256=digest)
        row["quality_eligible_now"] = row["eligible"] == "True" and row["snapshot_status"] == "valid" and status == "snapshot_content_verified"
    by_sample = defaultdict(list)
    for row in matched:
        key = row["biosample"] if known(row["biosample"]) else "unknown:" + row["genome_id"]
        by_sample[key].append(row)
    cells, drug_rows = [], []
    for drug in drugs:
        raw_counts, valid_counts, direct_counts, split_counts = Counter(), Counter(), Counter(), Counter()
        clone_classes, source_classes = defaultdict(set), defaultdict(set)
        raw_labeled_assemblies = sum(calls[r["genome_id"]][drug] in ({"0"}, {"1"}) for r in matched)
        conflicts = invalid = unknown = missing = quality_lost = 0
        for sample, records in sorted(by_sample.items()):
            states = set().union(*(calls[r["genome_id"]][drug] for r in records))
            status, label = "missing", ""
            labeled = [r for r in records if r["quality_eligible_now"] and calls[r["genome_id"]][drug] in ({"0"}, {"1"})]
            if not states:
                missing += 1
            elif {"0", "1"} <= states:
                status = "conflict"
                conflicts += 1
            elif states - {"0", "1"}:
                status = "invalid_nonbinary"
                invalid += 1
            elif sample.startswith("unknown:"):
                status = "unknown_biosample"
                unknown += 1
            else:
                label = "S" if states == {"0"} else "R"
                raw_counts[label] += 1
                if not labeled:
                    status = "quality_excluded"
                    quality_lost += 1
                else:
                    status = "available_not_selected"
                    valid_counts[label] += 1
                    if any(v[drug] == ("0" if label == "S" else "1") for r in labeled for v in by_accession.get(r["file_accession"], ())):
                        direct_counts[label] += 1
                    partitions = {r["inherited_split"] for r in labeled}
                    if len(partitions) != 1:
                        raise ValueError(f"BioSample spans inherited splits: {sample}")
                    partition = next(iter(partitions))
                    split_counts[f"{partition}_{label}"] += 1
                    for r in labeled:
                        clone_classes[r["near_clone_group"]].add(label)
                        if known(r["source_group"]):
                            source_classes[r["source_group"]].add(label)
            cells.append({"biosample": "" if sample.startswith("unknown:") else sample, "drug": drug, "status": status,
                          "label": label, "observed_states": json.dumps(sorted(states)),
                          "source_genome_ids": json.dumps([r["genome_id"] for r in records if calls[r["genome_id"]][drug]]),
                          "quality_eligible_labeled_genome_ids": json.dumps([r["genome_id"] for r in labeled]),
                          "inherited_split": labeled[0]["inherited_split"] if status == "available_not_selected" else ""})
        n_s, n_r = valid_counts["S"], valid_counts["R"]
        clone_support = {c: sum(c in v for v in clone_classes.values()) for c in ("S", "R")}
        source_support = {c: sum(c in v for v in source_classes.values()) for c in ("S", "R")}
        gate = n_s + n_r >= 500 and min(n_s, n_r) >= 50
        split_coverage = all(split_counts[f"{s}_{c}"] > 0 for s in SPLITS for c in ("S", "R"))
        necessary_clone = min(clone_support.values()) >= 3
        necessary_source = min(source_support.values()) >= 3
        if not necessary_clone or not necessary_source:
            decision = "not_currently_evaluable_three_way"
        elif gate:
            decision = "count_gate_passes_protocol_pending"
        else:
            decision = "exploratory_only_protocol_pending"
        drug_rows.append({"drug": drug, "raw_labeled_assemblies": raw_labeled_assemblies,
                          "raw_nonconflicting_biosamples": sum(raw_counts.values()), "raw_S": raw_counts["S"], "raw_R": raw_counts["R"],
                          "quality_eligible_biosamples": n_s + n_r, "S": n_s, "R": n_r,
                          "direct_file_biosamples_S": direct_counts["S"], "direct_file_biosamples_R": direct_counts["R"],
                          "alias_only_biosamples_S": n_s - direct_counts["S"], "alias_only_biosamples_R": n_r - direct_counts["R"],
                          "conflicting_biosamples": conflicts, "invalid_label_units": invalid, "unknown_biosample_units": unknown,
                          "missing_label_units": missing, "quality_excluded_biosamples": quality_lost,
                          "near_clone_groups": len(clone_classes), "near_clone_groups_S": clone_support["S"], "near_clone_groups_R": clone_support["R"],
                          "near_clone_groups_mixed_SR": sum(v == {"S", "R"} for v in clone_classes.values()),
                          "source_components": len(source_classes), "source_components_S": source_support["S"], "source_components_R": source_support["R"],
                          **{f"inherited_{s}_{c}": split_counts[f"{s}_{c}"] for s in SPLITS for c in ("S", "R")},
                          "passes_500_50_count_gate": gate, "inherited_splits_both_classes": split_coverage,
                          "necessary_clone_group_support": necessary_clone, "necessary_source_group_support": necessary_source,
                          "protocol_blocker": "insufficient_near_clone_class_groups" if not necessary_clone else "insufficient_source_class_groups" if not necessary_source else "count_gate_not_met" if not gate else "mapping_and_evaluation_protocol_pending",
                          "protocol_decision": decision})
    drug_rows.sort(key=lambda r: (-r["quality_eligible_biosamples"], r["drug"]))
    outputs = {"record_disposition.tsv": matched, "sample_drug_evidence.tsv": cells, "drug_feasibility.tsv": drug_rows}
    for name, records in outputs.items():
        write_tsv(output / name, records, list(records[0]))
    if input_hashes != {str(p): sha256(p) for p in inputs}:
        raise ValueError("Audit inputs changed during run")
    summary = {
        "status": "amr_feasibility_audit_complete", "owner": "gpt", "git_commit": contract["git_commit"],
        "matched_assemblies": len(matched), "label_rows": len(label_rows),
        "unique_known_biosamples": sum(not s.startswith("unknown:") for s in by_sample),
        "assemblies_with_unknown_biosample": sum(not known(r["biosample"]) for r in matched),
        "biosamples_with_multiple_assemblies": sum(len(v) > 1 for s, v in by_sample.items() if not s.startswith("unknown:")),
        "dispositions": dict(Counter(r["disposition"] for r in matched)),
        "legacy_loss_mechanisms": dict(Counter(r["loss_mechanism"] for r in matched)),
        "quality_exclusion_reasons_nonexclusive": dict(Counter(reason for r in matched for reason in r["exclusion_reason"].split(";") if reason)),
        "current_file_status": dict(Counter(r["current_file_status"] for r in matched)),
        "quality_eligible_assemblies_now": sum(r["quality_eligible_now"] for r in matched),
        "quality_eligible_near_clone_groups": len({r["near_clone_group"] for r in matched if r["quality_eligible_now"]}),
        "quality_eligible_inherited_splits": dict(Counter(r["inherited_split"] for r in matched if r["quality_eligible_now"])),
        "version_only_label_accessions": version_only,
        "drugs_passing_count_gate": [r["drug"] for r in drug_rows if r["passes_500_50_count_gate"]],
        "protocol_decisions": dict(Counter(r["protocol_decision"] for r in drug_rows)),
        "n_sample_drug_conflicts_observed": sum(r["conflicting_biosamples"] for r in drug_rows),
        "training_performed": False, "new_split_created": False,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "limitations": [
            "Counts describe versioned published binary labels, not newly verified AST measurements or epidemiological independence.",
            "Conflicts within the collapsed source CSV cells and experimental conditions are unobservable; zero observed conflict does not prove no AST conflict.",
            "Exact paired accession is database linkage, not independent proof of sequence identity; match fields are retained per assembly.",
            "Quality reuses frozen cohort rules and verified content, without new contamination or completeness estimation.",
            "Inherited splits use original full membership source/clone components, never newly assign samples; class coverage is a minimum condition, not statistical power.",
            "Neither a count gate nor necessary group support proves a valid protocol, model utility, cross-strain generalisation or protein generation.",
        ],
    }
    atomic_json(output / "summary.json", summary)
    atomic_json(output / "completion.json", {name: sha256(output / name) for name in ("summary.json", "contract.json", *outputs)})
    return summary
