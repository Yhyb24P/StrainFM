"""Versioned, resumable public assembly metadata for a species inventory."""

import gzip
import json
import re
import time
import urllib.parse
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .audit import sha256, write_tsv
from .metadata import known
from .snapshot import atomic_json, prepare_output


def flatten(report):
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


def fetch_batch(item):
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
    accessions = sorted({r["file_accession"] for r in inventory})
    if not all(re.fullmatch(r"GC[AF]_\d{9}\.\d+", a) for a in accessions):
        raise ValueError("Only public versioned assembly accession IDs are queried")
    contract = {"schema": 1, "inventory_sha256": sha256(inventory_path), "code_sha256": sha256(__file__), "accessions": len(accessions), "batch_size": 100}
    prepare_output(output, contract)
    requested = set(accessions)
    batches = [(i // 100, accessions[i:i + 100], output) for i in range(0, len(accessions), 100)]
    reports, requests = {}, []
    with ThreadPoolExecutor(max_workers=2) as pool:
        for i, (payload, request) in enumerate(pool.map(fetch_batch, batches), 1):
            requests.append(request)
            for report in payload.get("reports", []):
                if report["accession"] in requested:
                    reports[report["accession"]] = flatten(report)
            if i % 10 == 0 or i == len(batches):
                atomic_json(output / "progress.json", {"batches_completed": i, "batches_total": len(batches), "reports_found": len(reports), "updated_utc": datetime.now(timezone.utc).isoformat()})
                print(f"metadata batches {i}/{len(batches)}, exact accessions {len(reports)}", flush=True)
    columns = list(flatten({"accession": ""})) + ["report_status"]
    records = [{**reports[a], "report_status": "found"} if a in reports else {"file_accession": a, "report_status": "unknown"} for a in accessions]
    write_tsv(output / "assembly_metadata.tsv", records, columns)
    atomic_json(output / "requests.json", requests)
    coverage = {field: sum(known(r.get(field, "")) for r in records) for field in ("biosample", "bioproject", "collection_date", "geo_loc_name", "isolation_source", "contig_n50", "assembly_level")}
    result = {"status": "metadata_complete" if len(reports) == len(accessions) else "metadata_partial", "requested": len(accessions), "found_exact_version": len(reports), "unknown": len(accessions) - len(reports), "failed_batches": sum("error" in r for r in requests), "coverage_nonmissing": coverage, "assembly_levels": dict(Counter(r.get("assembly_level", "unknown") for r in records)), "metadata_sha256": sha256(output / "assembly_metadata.tsv"), "finished_utc": datetime.now(timezone.utc).isoformat(), "limitations": ["Metadata are fetched at the recorded times; assembly taxonomy and assay/source assertions are not independently established.", "Only exact requested accession versions are used; current or paired reports do not silently replace missing versions.", "Public BioSample attributes are saved locally; not automatically approved for public redistribution."]}
    atomic_json(output / "summary.json", result)
    return result
