"""Construct a species cohort from file evidence and assembly metadata."""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from .audit import sha256, write_tsv
from .metadata import known
from .snapshot import atomic_json


class UnionFind:
    def __init__(self, keys):
        self.parent = {k: k for k in keys}

    def find(self, key):
        while self.parent[key] != key:
            self.parent[key] = self.parent[self.parent[key]]
            key = self.parent[key]
        return key

    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)
            return True
        return False

    def groups(self):
        groups = defaultdict(list)
        for key in sorted(self.parent):
            groups[self.find(key)].append(key)
        return groups


def read_table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rank_quality(row):
    level = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    return (level.get(row.get("assembly_level"), 4), -int(row.get("contig_n50") or 0), int(row.get("contig_count") or 10**9), row["file_accession"])


def run_cohort(snapshot_path, metadata_path, config_path, output):
    snapshot_path, metadata_path, config_path, output = map(Path, (snapshot_path, metadata_path, config_path, output))
    if output.exists() or "raw" in output.resolve().parts:
        raise ValueError("Use a new non-raw output directory")
    config = json.loads(config_path.read_text())
    samples, metadata = read_table(snapshot_path), {r["file_accession"]: r for r in read_table(metadata_path)}
    records = []
    for sample in samples:
        report = metadata.get(sample["file_accession"], {})
        record = {**report, **sample}
        reasons = []
        if sample["status"] != "valid":
            reasons.append("file_" + sample["status"])
        if report.get("report_status") != "found":
            reasons.append("metadata_exact_version_unknown")
        else:
            species = report.get("submitted_species", "")
            if (species and species != config["scientific_name"]) or not report["organism_name"].startswith(config["scientific_name"]):
                reasons.append("metadata_species_mismatch")
            if report.get("assembly_status", "").lower() == "suppressed":
                reasons.append("assembly_suppressed")
            if report.get("taxonomy_check_status", "").upper() in ("FAILED", "FAIL"):
                reasons.append("taxonomy_check_failed")
            for field, minimum in (("contig_count", 1), ("contig_n50", config["min_contig_n50"])):
                if not known(report.get(field, "")):
                    reasons.append(field + "_unknown")
                elif int(report[field]) < minimum:
                    reasons.append(field + "_below_minimum")
            if known(report.get("contig_count", "")) and int(report["contig_count"]) > config["max_contigs"]:
                reasons.append("contig_count_above_maximum")
        if sample["status"] == "valid":
            if not config["min_bases"] <= int(sample["total_bases"]) <= config["max_bases"]:
                reasons.append("length_outside_range")
            if float(sample["ambiguous_fraction"]) > config["max_ambiguous_fraction"]:
                reasons.append("ambiguous_fraction_above_maximum")
            if report.get("total_sequence_length") and str(report["total_sequence_length"]) != sample["total_bases"]:
                reasons.append("length_disagrees_with_accession_metadata")
        record.update(eligible=not reasons, exclusion_reason=";".join(reasons))
        records.append(record)
    eligible = {r["genome_id"]: r for r in records if r["eligible"]}
    union = UnionFind(eligible)
    for field in ("biosample", "contig_multiset_sha256", "file_accession"):
        first = {}
        for gid, row in eligible.items():
            value = row.get(field, "")
            if known(value):
                if value in first:
                    union.union(first[value], gid)
                else:
                    first[value] = gid
    groups = union.groups()
    chosen = {}
    for group, ids in groups.items():
        representative = min(ids, key=lambda gid: rank_quality(eligible[gid]))
        for gid in ids:
            chosen[gid] = (group, representative)
    for row in records:
        if row["genome_id"] in chosen:
            group, representative = chosen[row["genome_id"]]
            row.update(identity_group=group, identity_representative=representative, retained_identity_representative=row["genome_id"] == representative)
        else:
            row.update(identity_group="", identity_representative="", retained_identity_representative=False)
    output.mkdir(parents=True)
    columns = list(dict.fromkeys(key for row in records for key in row))
    write_tsv(output / "cohort_manifest.tsv", records, columns)
    representatives = [r for r in records if r["retained_identity_representative"]]
    write_tsv(output / "identity_representatives.tsv", representatives, columns)
    source_root = Path(config["genome_root"])
    (output / "mash_inputs.txt").write_text("".join(str(source_root / row["target_relpath"]) + "\n" for row in representatives))
    counts = Counter(reason for row in records for reason in row["exclusion_reason"].split(";") if reason)
    summary = {"status": "identity_cohort_ready_for_sequence_grouping", "config": config, "input_sha256": {str(p): sha256(p) for p in (snapshot_path, metadata_path, config_path)}, "code_sha256": sha256(__file__), "records": len(records), "eligible_files": len(eligible), "excluded_files": len(records) - len(eligible), "exclusion_reasons_nonexclusive": dict(counts), "identity_groups": len(groups), "identity_duplicate_excess": len(eligible) - len(groups), "biosample_unknown_eligible": sum(not known(r.get("biosample", "")) for r in eligible.values()), "bioproject_unknown_representatives": sum(not known(r.get("bioproject", "")) for r in representatives), "limitations": ["Identity groups merge eligible assemblies by exact FASTA-record digest or BioSample, not by genomic near-clone distance.", "NCBI assembly metadata support filtering; no new annotation or contamination estimation was run.", "Pretraining overlap remains unknown until a usable upstream membership list is obtained."]}
    atomic_json(output / "summary.json", summary)
    return summary
