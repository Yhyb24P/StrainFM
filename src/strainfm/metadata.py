"""NCBI assembly metadata: sampled audits and resumable inventory queries."""

import gzip
import hashlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .audit import rows, sha256, write_tsv
from .snapshot import atomic_json, prepare_output


def known(value):
    return str(value).strip().lower() not in {"", "none", "na", "n/a", "not applicable", "not collected", "not provided", "missing", "unknown", "not available"}


def run_metadata_audit(sample_path, output):
    sample_path, output = Path(sample_path), Path(output)
    if "raw" in output.resolve().parts:
        raise ValueError("Output must be outside raw directories")
    sample = list(rows(sample_path, ("genome_id", "file_accession", "species_id")))
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
            fields = _flatten(report)
            info = report.get("assembly_info", {})
            biosample = info.get("biosample", {})
            attrs = {a["name"]: a.get("value", "") for a in biosample.get("attributes", [])}
            date = fields["collection_date"]
            location = fields["geo_loc_name"] or attrs.get("geographic location", "")
            isolation = fields["isolation_source"]
            paired = fields["paired_accession"]
            record.update({"reported_organism": fields["organism_name"],
                           "tax_id": fields["tax_id"],
                           "assembly_level": fields["assembly_level"],
                           "assembly_status": fields["assembly_status"],
                           "paired_accession": paired,
                           "original_id_is_returned_or_paired": row["genome_id"] in (accession, paired),
                           "biosample": fields["biosample"],
                           "bioproject": fields["bioproject"],
                           "collection_date_known": known(date), "location_known": known(location),
                           "isolation_source_known": known(isolation), "sequencing_tech_known": known(fields["sequencing_tech"]),
                           "ncbi_contig_count": fields["contig_count"],
                           "ncbi_contig_n50": fields["contig_n50"],
                           "ncbi_scaffold_count": fields["scaffold_count"],
                           "ncbi_scaffold_n50": fields["scaffold_n50"],
                           "ncbi_total_length": fields["total_sequence_length"],
                           "fasta_total_length_matches": str(fields["total_sequence_length"]) == row.get("total_bases", "unavailable"),
                           "annotated_protein_coding_genes": fields["protein_coding_genes"]})
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


def _flatten(report):
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    bio = info.get("biosample", {})
    attrs = {a["name"]: a.get("value", "") for a in bio.get("attributes", [])}
    ani = report.get("average_nucleotide_identity", {})
    return {
        "file_accession": report["accession"],
        "current_accession": report.get("current_accession", ""),
        "paired_accession": report.get("paired_accession") or info.get("paired_assembly", {}).get("accession", ""),
        "organism_name": report.get("organism", {}).get("organism_name", ""),
        "tax_id": report.get("organism", {}).get("tax_id", ""),
        "submitted_species": ani.get("submitted_species", ""),
        "taxonomy_check_status": ani.get("taxonomy_check_status", ""),
        "assembly_level": info.get("assembly_level", ""),
        "assembly_status": info.get("assembly_status", ""),
        "biosample": bio.get("accession", ""),
        "bioproject": info.get("bioproject_accession", ""),
        "collection_date": bio.get("collection_date") or attrs.get("collection_date", ""),
        "geo_loc_name": bio.get("geo_loc_name") or attrs.get("geo_loc_name", ""),
        "isolation_source": bio.get("isolation_source") or attrs.get("isolation_source", ""),
        "sequencing_tech": info.get("sequencing_tech", ""),
        "total_sequence_length": stats.get("total_sequence_length", ""),
        "contig_count": stats.get("number_of_contigs", ""),
        "contig_n50": stats.get("contig_n50", ""),
        "scaffold_count": stats.get("number_of_scaffolds", ""),
        "scaffold_n50": stats.get("scaffold_n50", ""),
        "protein_coding_genes": report.get("annotation_info", {}).get("stats", {}).get("gene_counts", {}).get("protein_coding", ""),
    }


def _fetch_batch(item):
    index, accessions, output = item
    path = output / f"batch_{index:04d}.json.gz"
    metadata = output / f"batch_{index:04d}.request.json"
    if path.exists() and metadata.exists():
        request = json.loads(metadata.read_text())
        if request["requested_accessions"] != accessions or sha256(path) != request["response_gzip_sha256"]:
            raise ValueError("Cached metadata batch changed")
        return json.loads(gzip.decompress(path.read_bytes())), request
    base = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
    url = base + urllib.parse.quote(",".join(accessions), safe=",._") + "/dataset_report?filters.assembly_version=all_assemblies&page_size=1000"
    request = {"url": url, "requested_accessions": accessions}
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=45) as response:
                payload = json.load(response)
            if payload.get("next_page_token"):
                raise ValueError("Unexpected pagination for <=100 requested assemblies")
            if "reports" not in payload:
                raise ValueError("Response has no reports field")
            content = gzip.compress(json.dumps(payload, separators=(",", ":")).encode(), mtime=0)
            temporary = path.with_suffix(".tmp")
            temporary.write_bytes(content)
            temporary.replace(path)
            request.update(accessed_utc=datetime.now(timezone.utc).isoformat(), response_gzip_sha256=sha256(path))
            atomic_json(metadata, request)
            return payload, request
        except (OSError, ValueError) as exc:
            if attempt == 2:
                request.update(error=str(exc), failed_utc=datetime.now(timezone.utc).isoformat())
                return {"reports": []}, request
            time.sleep(2 ** attempt)


def run_bulk_metadata(inventory_path, output):
    inventory_path, output = Path(inventory_path), Path(output)
    if "raw" in output.resolve().parts:
        raise ValueError("Output must be outside raw")
    inventory = json.loads(inventory_path.read_text())
    if (not isinstance(inventory, list) or not inventory
            or any(not isinstance(record, dict) or not isinstance(record.get("file_accession"), str)
                   for record in inventory)):
        raise ValueError("Inventory must contain records with versioned file_accession strings")
    accessions = sorted({record["file_accession"] for record in inventory})
    if not all(re.fullmatch(r"GC[AF]_\d{9}\.\d+", a) for a in accessions):
        raise ValueError("Only public versioned assembly accession IDs are queried")
    contract = {"schema": 1, "inventory_sha256": sha256(inventory_path), "code_sha256": sha256(__file__), "accessions": len(accessions), "batch_size": 100}
    prepare_output(output, contract)
    requested = set(accessions)
    batches = [(i // 100, accessions[i:i + 100], output) for i in range(0, len(accessions), 100)]
    reports, requests = {}, []
    with ThreadPoolExecutor(max_workers=2) as pool:
        for i, (payload, request) in enumerate(pool.map(_fetch_batch, batches), 1):
            requests.append(request)
            for report in payload.get("reports", []):
                if report["accession"] in requested:
                    reports[report["accession"]] = _flatten(report)
            if i % 10 == 0 or i == len(batches):
                atomic_json(output / "progress.json", {"batches_completed": i, "batches_total": len(batches), "reports_found": len(reports), "updated_utc": datetime.now(timezone.utc).isoformat()})
                print(f"metadata batches {i}/{len(batches)}, exact accessions {len(reports)}", flush=True)
    columns = list(_flatten({"accession": ""})) + ["report_status"]
    records = [{**reports[a], "report_status": "found"} if a in reports else {"file_accession": a, "report_status": "unknown"} for a in accessions]
    write_tsv(output / "assembly_metadata.tsv", records, columns)
    atomic_json(output / "requests.json", requests)
    coverage = {field: sum(known(r.get(field, "")) for r in records) for field in ("biosample", "bioproject", "collection_date", "geo_loc_name", "isolation_source", "contig_n50", "assembly_level")}
    result = {"status": "metadata_complete" if len(reports) == len(accessions) else "metadata_partial", "requested": len(accessions), "found_exact_version": len(reports), "unknown": len(accessions) - len(reports), "failed_batches": sum("error" in r for r in requests), "coverage_nonmissing": coverage, "assembly_levels": dict(Counter(r.get("assembly_level", "unknown") for r in records)), "metadata_sha256": sha256(output / "assembly_metadata.tsv"), "finished_utc": datetime.now(timezone.utc).isoformat(), "limitations": ["Metadata are fetched at the recorded times; assembly taxonomy and assay/source assertions are not independently established.", "Only exact requested accession versions are used; current or paired reports do not silently replace missing versions.", "Public BioSample attributes are saved locally; not automatically approved for public redistribution."]}
    atomic_json(output / "summary.json", result)
    return result
