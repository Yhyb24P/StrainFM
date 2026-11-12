"""Link a fixed BV-BRC AST snapshot to samples and candidate genome files."""

import csv
import json
import lzma
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path

from .audit import safe_path, sha256, write_tsv
from .cohort import rank_quality, read_table
from .snapshot import atomic_json

DRUGS = ("ciprofloxacin", "gentamicin")
BIOSAMPLE = re.compile(r"SAM[NED][A-Z]?\d+")


def _read_pages(manifest, source, hashes):
    hashes[str(manifest)] = sha256(manifest)
    records = []
    for page in json.loads(manifest.read_text()):
        path = Path(page["path"])
        if not path.resolve().is_relative_to(source.resolve()):
            raise ValueError(f"Response outside source snapshot: {path}")
        digest = sha256(path)
        if digest != page["sha256"]:
            raise ValueError(f"Response hash mismatch: {path}")
        values = json.loads(path.read_text())
        if not isinstance(values, list) or len(values) != page["rows"]:
            raise ValueError(f"Response row count mismatch: {path}")
        hashes[str(path)] = digest
        records.extend(values)
    return records


def _aliases(source, hashes):
    aliases = {}
    for path in sorted((source / "ena_aliases").glob("batch_*.xml")):
        request = path.with_suffix(".request.json")
        hashes[str(request)] = sha256(request)
        digest = sha256(path)
        if digest != json.loads(request.read_text())["sha256"]:
            raise ValueError(f"ENA response hash mismatch: {path}")
        hashes[str(path)] = digest
        for node in ET.fromstring(path.read_bytes()).findall("SAMPLE"):
            samples = {x.text for x in node.findall("IDENTIFIERS/EXTERNAL_ID")
                       if x.get("namespace") == "BioSample"}
            if len(samples) != 1 or not BIOSAMPLE.fullmatch(next(iter(samples)) or ""):
                raise ValueError(f"Ambiguous ENA BioSample in {path}")
            sample = next(iter(samples))
            keys = {node.get("accession"), node.findtext("IDENTIFIERS/SECONDARY_ID")}
            for key in keys - {None, ""}:
                if key in aliases and aliases[key] != sample:
                    raise ValueError(f"Conflicting ENA alias: {key}")
                aliases[key] = sample
    return aliases


def phenotype_state(records):
    states = {r.get("resistant_phenotype", "") for r in records} - {"", None}
    if not states:
        return "unknown"
    if len(states) != 1:
        return "conflict"
    return {"Susceptible": "S", "Resistant": "R", "Intermediate": "I"}.get(states.pop(), "other")


def link_samples(records, metadata, aliases, cohort):
    """Never transfer a label across samples or interpret an unreported AST call."""
    by_sample, by_accession = defaultdict(dict), defaultdict(dict)
    for row in cohort:
        if BIOSAMPLE.fullmatch(row.get("biosample", "")):
            by_sample[row["biosample"]][row["genome_id"]] = row
        for field in ("file_accession", "current_accession", "paired_accession"):
            if row.get(field):
                by_accession[row[field]][row["genome_id"]] = row
    genome_links = {}
    for gid, row in metadata.items():
        raw = row.get("biosample_accession") or ""
        sample = aliases.get(raw, raw)
        sample = sample if BIOSAMPLE.fullmatch(sample) else ""
        accession_hits = by_accession.get(row.get("assembly_accession", ""), {})
        conflict = bool(sample and any(r.get("biosample") and r["biosample"] != sample
                                      for r in accession_hits.values()))
        genome_links[gid] = (sample, conflict)
    grouped = defaultdict(list)
    for row in records:
        sample, _ = genome_links[row["genome_id"]]
        grouped[(sample or "unresolved:" + row["genome_id"], row["antibiotic"])].append(row)
    units = []
    for (key, drug), values in sorted(grouped.items()):
        gids = sorted({r["genome_id"] for r in values})
        sample = "" if key.startswith("unresolved:") else key
        units.append({
            "sample_key": key, "biosample": sample, "drug": drug,
            "state": phenotype_state(values), "record_count": len(values),
            "unknown_record_count": sum(not r.get("resistant_phenotype") for r in values),
            "bvbrc_genome_ids": ";".join(gids),
            "record_ids": ";".join(sorted(r["id"] for r in values)),
            "identity_conflict": any(genome_links[g][1] for g in gids),
        })
    samples = []
    for sample in sorted({r["biosample"] for r in units if r["biosample"] and r["state"] in ("S", "R")}):
        sample_units = [r for r in units if r["biosample"] == sample]
        conflict = any(r["identity_conflict"] for r in sample_units)
        eligible = [r for r in by_sample.get(sample, {}).values()
                    if r.get("eligible") == "True" and r.get("status") == "valid"]
        chosen = min(eligible, key=rank_quality) if eligible and not conflict else None
        row = {"genome_id": sample, "biosample": sample,
               "action": "quarantine_identity_conflict" if conflict else "reuse_local" if chosen else "download_atb",
               "local_genome_id": chosen["genome_id"] if chosen else "",
               "file_accession": chosen["file_accession"] if chosen else "",
               "local_target_relpath": chosen["target_relpath"] if chosen else "",
               "expected_local_sha256": chosen["file_sha256"] if chosen else "",
               "source_bvbrc_genome_ids": ";".join(sorted({g for r in sample_units for g in r["bvbrc_genome_ids"].split(";")}))}
        for drug in DRUGS:
            observed = [r["state"] for r in sample_units if r["drug"] == drug]
            row[drug] = observed[0] if observed else "unknown"
        samples.append(row)
    return units, samples


def apply_atb_metadata(samples, metadata_path):
    """Gate candidate downloads on the fixed ATB assembly table, not URL existence."""
    selected = {r["biosample"]: r for r in samples if r["action"] == "download_atb"}
    found = set()
    opener = lzma.open if metadata_path.suffix == ".xz" else open
    with opener(metadata_path, "rt") as handle:
        # ATB explicitly requires literal quote parsing for this exported table.
        reader = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
        required = {"sample_accession", "run_accession", "assembly_accession", "dataset", "asm_pipe_filter",
                    "asm_fasta_on_osf", "sylph_species", "sylph_filter", "hq_filter", "aws_url"}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Missing ATB assembly metadata columns")
        for meta in reader:
            key = meta["sample_accession"]
            if key not in selected:
                continue
            if key in found:
                raise ValueError(f"Duplicate ATB sample: {key}")
            found.add(key)
            row = selected[key]
            row.update(atb_dataset=meta["dataset"], atb_run_accessions=meta["run_accession"],
                       atb_assembly_accession=meta["assembly_accession"], atb_species=meta["sylph_species"])
            reasons = [f"{field}={meta[field]}" for field in ("asm_pipe_filter", "sylph_filter", "hq_filter")
                       if meta[field] != "PASS"]
            if meta["asm_fasta_on_osf"] != "1":
                reasons.append("assembly_unavailable")
            if meta["sylph_species"] != "Escherichia coli":
                reasons.append("species_not_Escherichia_coli")
            if not meta["run_accession"] or meta["run_accession"] == "NA":
                reasons.append("assembly_run_unknown")
            if meta["aws_url"] != "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/" + key + ".fa.gz":
                reasons.append("unexpected_assembly_url")
            row["atb_exclusion_reason"] = ";".join(reasons)
            if reasons:
                row["action"] = "atb_metadata_excluded"
    for key in selected.keys() - found:
        selected[key].update(action="atb_metadata_missing", atb_exclusion_reason="sample_not_in_metadata")
    for row in samples:
        for key in ("atb_dataset", "atb_run_accessions", "atb_assembly_accession", "atb_species", "atb_exclusion_reason"):
            row.setdefault(key, "")


def apply_ena_runs(samples, source, hashes):
    """Retain all study links and quarantine conflicting source taxonomy."""
    manifest = source / "manifest.json"
    hashes[str(manifest)] = sha256(manifest)
    queried, seen_runs = set(), set()
    by_sample = defaultdict(list)
    for page in json.loads(manifest.read_text()):
        path = Path(page["path"])
        if not path.resolve().is_relative_to(source.resolve()):
            raise ValueError(f"ENA run response outside snapshot: {path}")
        if sha256(path) != page["sha256"]:
            raise ValueError(f"ENA run response hash mismatch: {path}")
        hashes[str(path)] = page["sha256"]
        requested = set(page["requested_samples"])
        if queried & requested:
            raise ValueError("Duplicate ENA sample query")
        queried.update(requested)
        records = read_table(path)
        if len(records) != page["returned_rows"]:
            raise ValueError(f"ENA run row count mismatch: {path}")
        for record in records:
            sample, run = record["sample_accession"], record["run_accession"]
            if sample not in requested or not run or run in seen_runs:
                raise ValueError(f"Unexpected or duplicate ENA run: {run}")
            seen_runs.add(run)
            by_sample[sample].append(record)
    for row in samples:
        sample = row["biosample"]
        records = by_sample[sample]
        names = {r["scientific_name"] for r in records} - {""}
        conflict = any(n != "Escherichia coli" and not n.startswith("Escherichia coli ") for n in names)
        row.update(ena_run_accessions=";".join(sorted({r["run_accession"] for r in records})),
                   ena_study_accessions=";".join(sorted({r["study_accession"] for r in records} - {""})),
                   ena_species=";".join(sorted(names)),
                   ena_query_status="records" if records else "no_records" if sample in queried else "not_queried",
                   ena_taxonomy_conflict=conflict)
        if conflict and row["action"] in ("download_atb", "reuse_local"):
            row["action"] = "quarantine_source_taxonomy_conflict"


def run_amr_sources(source, cohort_path, local_genome_root, output, *, atb_assemblies=None, ena_runs=None):
    source, cohort_path, local_genome_root, output = map(Path, (source, cohort_path, local_genome_root, output))
    if output.exists() or "raw" in output.resolve().parts or output.resolve().is_relative_to(source.resolve()):
        raise ValueError("Use a new output directory outside raw and the source snapshot")
    hashes = {str(cohort_path): sha256(cohort_path)}
    records = []
    for drug in DRUGS:
        manifest = source / "bvbrc" / f"{drug}.manifest.json"
        values = _read_pages(manifest, source, hashes)
        expected = {int(p["content_range"].split("/")[-1]) for p in json.loads(manifest.read_text())}
        if expected != {len(values)}:
            raise ValueError(f"Incomplete AST snapshot: {drug}")
        for row in values:
            if (row.get("taxon_id") != 562 or row.get("antibiotic") != drug
                    or row.get("evidence") != "Laboratory Method" or row.get("public") is not True
                    or row.get("computational_method")):
                raise ValueError(f"Record outside measured-AST contract: {row.get('id')}")
        records.extend(values)
    if len({r["id"] for r in records}) != len(records):
        raise ValueError("Duplicate AST record IDs")
    meta_rows = _read_pages(source / "bvbrc" / "metadata" / "manifest.json", source, hashes)
    metadata = {r["genome_id"]: r for r in meta_rows}
    if len(metadata) != len(meta_rows) or {r["genome_id"] for r in records} != set(metadata):
        raise ValueError("Metadata must cover each AST genome exactly once")
    units, samples = link_samples(records, metadata, _aliases(source, hashes), read_table(cohort_path))
    if not samples:
        raise ValueError("No samples with an interpretable binary AST call")
    if atb_assemblies is not None:
        atb_assemblies = Path(atb_assemblies)
        hashes[str(atb_assemblies)] = sha256(atb_assemblies)
        apply_atb_metadata(samples, atb_assemblies)
    if ena_runs is not None:
        apply_ena_runs(samples, Path(ena_runs), hashes)
    downloads = []
    for row in samples:
        row["local_fasta_path"] = ""
        if row["action"] == "reuse_local":
            path = safe_path(local_genome_root, row["local_target_relpath"])
            if not path.is_file() or sha256(path) != row["expected_local_sha256"]:
                row["action"] = "quarantine_local_file_changed"
            else:
                row["local_fasta_path"] = str(path)
        if row["action"] == "download_atb":
            downloads.append({"genome_id": row["genome_id"], "biosample": row["biosample"],
                              "target_relpath": row["biosample"] + ".fa.gz",
                              "url": "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/" + row["biosample"] + ".fa.gz"})
    output.mkdir(parents=True)
    write_tsv(output / "sample_drug_evidence.tsv", units, list(units[0]))
    write_tsv(output / "sample_genomes.tsv", samples, list(samples[0]))
    write_tsv(output / "download_manifest.tsv", downloads, ["genome_id", "biosample", "target_relpath", "url"])
    summary = {"status": "linked_candidates_not_model_ready", "input_sha256": hashes,
               "code_sha256": sha256(__file__), "local_genome_root": str(local_genome_root),
               "atb_metadata_checked": atb_assemblies is not None,
               "ena_runs_checked": ena_runs is not None,
               "ast_records": len(records), "binary_candidate_samples": len(samples),
               "actions": dict(Counter(r["action"] for r in samples)),
               "drug_states": {d: dict(Counter(r["state"] for r in units if r["drug"] == d and r["biosample"])) for d in DRUGS},
               "limits": ["Source scope: public Laboratory Method records, taxon_id 562, two drugs.",
                          "Unresolved samples, I, unknown calls and label conflicts stay in the evidence table.",
                          "ATB URLs are candidates; release metadata, file validation, QC and grouping remain required.",
                          "No cross-sample label transfer or version-prefix accession match; no data split or model processing."]}
    atomic_json(output / "summary.json", summary)
    return summary
