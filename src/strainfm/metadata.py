"""Small, sequential NCBI metadata audit; no sequence downloads."""

import csv
import hashlib
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv


def known(value):
    return str(value).strip().lower() not in {"", "none", "na", "n/a", "not applicable", "not collected", "not provided", "missing", "unknown", "not available"}


def run_metadata_audit(sample_path, output):
    sample_path, output = Path(sample_path), Path(output)
    if "raw" in output.resolve().parts:
        raise ValueError("Output must be outside raw directories")
    with sample_path.open(newline="") as handle:
        sample = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted({row["file_accession"] for row in sample})
    if not accessions or len(accessions) > 200:
        raise ValueError("Metadata audit is limited to 1–200 sampled accessions")
    if not all(re.fullmatch(r"GC[AF]_\d{9}\.\d+", a) for a in accessions):
        raise ValueError("Only public versioned assembly accessions may be queried")
    input_digest = sha256(sample_path)
    output.mkdir(parents=True, exist_ok=False)
    reports, requests, failures = {}, [], []
    for start in range(0, len(accessions), 20):
        batch = accessions[start:start + 20]
        url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/" + urllib.parse.quote(",".join(batch), safe=",._") + "/dataset_report?page_size=100"
        item = {"url": url, "requested_accessions": batch, "accessed_utc": datetime.now(timezone.utc).isoformat()}
        try:
            with urllib.request.urlopen(url, timeout=25) as response:
                content = response.read()
            payload = json.loads(content)
            name = f"ncbi_batch_{start // 20:02d}.json"
            (output / name).write_bytes(content)
            item.update(file=name, sha256=hashlib.sha256(content).hexdigest())
            if payload.get("next_page_token"):
                item["warning"] = "Unexpected pagination; audit is incomplete"
                failures.append(item["warning"])
            for report in payload.get("reports", []):
                reports[report["accession"]] = report
        except (OSError, ValueError, urllib.error.URLError) as exc:
            item["error"] = str(exc)
            failures.append(str(exc))
        requests.append(item)
    records = []
    for row in sample:
        accession = row["file_accession"]
        report = reports.get(accession)
        record = {"genome_id": row["genome_id"], "file_accession": accession, "species_from_manifest": row["species_id"], "report_found": report is not None}
        if report is not None:
            info, stats = report.get("assembly_info", {}), report.get("assembly_stats", {})
            biosample = info.get("biosample", {})
            attrs = {a["name"]: a.get("value", "") for a in biosample.get("attributes", [])}
            date = biosample.get("collection_date") or attrs.get("collection_date", "")
            location = biosample.get("geo_loc_name") or attrs.get("geo_loc_name") or attrs.get("geographic location", "")
            isolation = biosample.get("isolation_source") or attrs.get("isolation_source", "")
            paired = report.get("paired_accession") or info.get("paired_assembly", {}).get("accession", "")
            record.update({"reported_organism": report.get("organism", {}).get("organism_name", ""),
                           "tax_id": report.get("organism", {}).get("tax_id", ""),
                           "assembly_level": info.get("assembly_level", ""),
                           "assembly_status": info.get("assembly_status", ""),
                           "paired_accession": paired,
                           "original_id_is_returned_or_paired": row["genome_id"] in (accession, paired),
                           "biosample": biosample.get("accession", ""),
                           "bioproject": info.get("bioproject_accession", ""),
                           "collection_date_known": known(date), "location_known": known(location),
                           "isolation_source_known": known(isolation), "sequencing_tech_known": known(info.get("sequencing_tech", "")),
                           "ncbi_contig_count": stats.get("number_of_contigs", ""),
                           "ncbi_contig_n50": stats.get("contig_n50", ""),
                           "ncbi_scaffold_count": stats.get("number_of_scaffolds", ""),
                           "ncbi_scaffold_n50": stats.get("scaffold_n50", ""),
                           "ncbi_total_length": stats.get("total_sequence_length", ""),
                           "fasta_total_length_matches": str(stats.get("total_sequence_length", "")) == row.get("total_bases", "unavailable"),
                           "annotated_protein_coding_genes": report.get("annotation_info", {}).get("stats", {}).get("gene_counts", {}).get("protein_coding", "")})
        records.append(record)
    columns = list(dict.fromkeys(key for record in records for key in record))
    write_tsv(output / "sample_metadata.tsv", records, columns)
    (output / "requests.json").write_text(json.dumps(requests, indent=2) + "\n")
    by_species = []
    for sp in sorted({r["species_from_manifest"] for r in records}):
        group = [r for r in records if r["species_from_manifest"] == sp]
        found = [r for r in group if r["report_found"]]
        summary = {"species": sp, "requested": len(group), "found": len(found), "assembly_levels": dict(Counter(r["assembly_level"] for r in found)),
                   "known_bioprojects": len({r["bioproject"] for r in found if known(r["bioproject"])}),
                   "missing_bioproject_records": sum(not known(r["bioproject"]) for r in found)}
        for key in ("collection_date_known", "location_known", "isolation_source_known", "sequencing_tech_known", "fasta_total_length_matches", "original_id_is_returned_or_paired"):
            summary[key] = sum(r[key] for r in found)
        summary["annotated_protein_coding_genes_gt6000"] = sum(int(r["annotated_protein_coding_genes"] or 0) > 6000 for r in found)
        summary["missing_protein_count"] = sum(r["annotated_protein_coding_genes"] == "" for r in found)
        by_species.append(summary)
    summary = {"owner": "Codex current session", "finished_utc": datetime.now(timezone.utc).isoformat(), "argv": sys.argv,
               "input_sha256": {str(sample_path): input_digest}, "code_sha256": {str(p): sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
               "status": "completed" if not failures and set(accessions) <= reports.keys() else "partial",
               "failures": failures, "species": by_species,
               "limitations": ["Metadata snapshot is not sequence-identity or upstream checksum verification.", "Paired accession validates a database relationship, not identical sequences.", "Missingness uses documented flattened fields and standard attribute names; unusual names can be missed.", "BioProject count does not establish epidemiological independence.", "Protein counts refer to deposited annotations, not a unified Bakta run or actually retained model tokens.", "FASTA records may be scaffolds; record N50 must not be interpreted as contig N50."]}
    if sha256(sample_path) != input_digest:
        summary["status"] = "input_changed"
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
