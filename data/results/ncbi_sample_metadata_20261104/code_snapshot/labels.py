"""Label availability, not phenotype prediction or evaluation cohort selection."""

import csv
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .audit import rows, sha256, write_tsv


def run_label_audit(labels_path, manifest_path, output):
    labels_path, manifest_path, output = map(Path, (labels_path, manifest_path, output))
    if "raw" in output.resolve().parts:
        raise ValueError("Output must be outside raw directories")
    if output.exists():
        raise FileExistsError(output)
    inputs = {str(p): sha256(p) for p in (labels_path, manifest_path)}
    original_ids, file_ids = defaultdict(set), defaultdict(set)
    for row in rows(manifest_path, ("genome_id", "file_accession", "species_id", "kingdom")):
        if row["kingdom"] == "bacteria":
            original_ids[row["genome_id"]].add(row["species_id"])
            file_ids[row["file_accession"]].add(row["species_id"])
    seen = set()
    global_counts, match_counts = Counter(), Counter()
    species_counts = defaultdict(Counter)
    by_drug = defaultdict(Counter)
    matched = []
    with labels_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "genome_name" not in (reader.fieldnames or []):
            raise ValueError("Missing genome_name")
        drugs = [k for k in reader.fieldnames if k != "genome_name"]
        for row in reader:
            gid = row["genome_name"]
            if not gid or gid in seen or None in row:
                raise ValueError(f"Malformed or duplicate label row: {gid}")
            seen.add(gid)
            values = {}
            for drug in drugs:
                value = row[drug]
                if value is None:
                    raise ValueError(f"Missing CSV cell: {gid}, {drug}")
                if value == "":
                    continue
                if value not in ("0", "0.0", "1", "1.0"):
                    raise ValueError(f"Unexpected binary label: {value}")
                values[drug] = int(float(value))
            global_counts["label_rows"] += 1
            global_counts["rows_with_any_label"] += bool(values)
            global_counts["observed_cells"] += len(values)
            global_counts["positive_cells"] += sum(values.values())
            categories = original_ids.get(gid, set()) | file_ids.get(gid, set())
            if not categories:
                match_counts["unmatched"] += 1
                continue
            if len(categories) > 1:
                match_counts["conflicting_species"] += 1
                continue
            sp = next(iter(categories))
            match_type = "both" if gid in original_ids and gid in file_ids else ("genome_id_only" if gid in original_ids else "file_accession_only")
            match_counts[match_type] += 1
            species_counts[sp]["matched_accessions"] += 1
            species_counts[sp]["accessions_with_any_label"] += bool(values)
            species_counts[sp]["observed_cells"] += len(values)
            matched.append({"label_accession": gid, "species_from_manifest": sp, "match_type": match_type, "observed_labels": len(values)})
            for drug, value in values.items():
                by_drug[(sp, drug)][value] += 1
    per_species = [{"species": sp, **dict(count)} for sp, count in sorted(species_counts.items())]
    per_drug = [{"species": sp, "drug": drug, "n_observed": count[0] + count[1], "n_label_0": count[0], "n_label_1": count[1], "positive_fraction": count[1] / (count[0] + count[1])} for (sp, drug), count in sorted(by_drug.items())]
    for path, digest in inputs.items():
        if sha256(path) != digest:
            raise ValueError(f"Input changed during audit: {path}")
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "species_label_coverage.tsv", per_species, ["species", "matched_accessions", "accessions_with_any_label", "observed_cells"])
    write_tsv(output / "species_drug_counts.tsv", per_drug, ["species", "drug", "n_observed", "n_label_0", "n_label_1", "positive_fraction"])
    write_tsv(output / "matched_accessions.tsv", matched, ["label_accession", "species_from_manifest", "match_type", "observed_labels"])
    summary = {"status": "label_availability_only_not_independent_evaluation", "owner": "Codex current session", "finished_utc": datetime.now(timezone.utc).isoformat(), "argv": sys.argv, "input_sha256": inputs,
               "code_sha256": {str(p): sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
               "global": {**dict(global_counts), "drug_columns": len(drugs), "missing_cells": len(seen) * len(drugs) - global_counts["observed_cells"]}, "matching": dict(match_counts),
               "limitations": ["Exact versioned accession match only; no GCF/GCA rewriting or version stripping.", "Species taken from local download manifest, not independently verified taxonomy.", "File presence and paired-assembly identity not verified by this command.", "Label source, assay method, repeated isolate independence and pretraining overlap need separate audit.", "Multiple drug labels on one genome are not independent genome observations."]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
