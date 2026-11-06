"""Read-only audits. No downloads, sequence annotation or quality-tier assignment."""

import csv
import gzip
import hashlib
import json
import os
import platform
import sys
import zlib
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def rows(path, required):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError(f"Missing required columns in {path}: {required}")
        for number, row in enumerate(reader, 2):
            extra = row.pop(None, [])
            if any(extra) or any(row.get(k) in (None, "") for k in required):
                raise ValueError(f"Malformed row {number} in {path}")
            row["_trailing_empty_fields"] = len(extra)
            yield row


def safe_path(root, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    target = (root / path).resolve()
    if not target.is_relative_to(root.resolve()):
        raise ValueError(f"Path escapes genome root: {relative}")
    return target


def sample_rank(seed, genome_id):
    return hashlib.sha256(f"{seed}\t{genome_id}".encode()).hexdigest()


def fasta_stats(path):
    """Fully decode gzip; N50 denominator is all bases, including ambiguous bases.

    Sequence digest ignores case, wrapping, headers and contig order, but not
    reverse complements, circular rotations or different contig boundaries.
    """
    path = Path(path)
    before = path.stat()
    compressed = path.read_bytes()
    compressed_digest = hashlib.sha256(compressed).hexdigest()
    try:
        content = gzip.decompress(compressed) if path.suffix == ".gz" else compressed
    except zlib.error as exc:
        raise ValueError(f"Invalid gzip data: {exc}") from exc
    lengths, contig_digests = [], []
    counts = Counter()
    if not content.strip() or not content.lstrip().startswith(b">"):
        raise ValueError("Empty FASTA or sequence before first header")
    allowed = b"ACGTRYSWKMBDHVN"
    for record in content.lstrip()[1:].split(b"\n>"):
        lines = record.splitlines()
        sequence = b"".join(line.strip() for line in lines[1:]).upper()
        if not sequence:
            raise ValueError("Empty FASTA record")
        if sequence.translate(None, allowed):
            raise ValueError("Non-IUPAC DNA characters")
        lengths.append(len(sequence))
        contig_digests.append(hashlib.sha256(sequence).hexdigest())
        counts.update({chr(base): sequence.count(base) for base in allowed})
    total = sum(lengths)
    cumulative = 0
    for contig_length in sorted(lengths, reverse=True):
        cumulative += contig_length
        if cumulative * 2 >= total:
            n50 = contig_length
            break
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Input changed during audit")
    acgt = sum(counts[base] for base in "ACGT")
    return {
        "bytes": before.st_size,
        "file_sha256": compressed_digest,
        "contig_multiset_sha256": hashlib.sha256("\n".join(sorted(contig_digests)).encode()).hexdigest(),
        "sequence_count": len(lengths),
        "total_bases": total,
        "n50_bases": n50,
        "max_contig_bases": max(lengths),
        "n_fraction": counts["N"] / total,
        "ambiguous_fraction": (total - acgt) / total,
        "gc_fraction_acgt": (counts["G"] + counts["C"]) / acgt if acgt else None,
    }


def write_tsv(path, records, columns):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def run_audit(config_path, output):
    config_path, output = Path(config_path), Path(output)
    config = json.loads(config_path.read_text())
    original = Path(config["manifest"])
    mapping = Path(config["download_manifest"])
    root = Path(config["genome_root"])
    species = config["species"]
    if len(species) != len(set(species)) or config["sample_per_species"] < 1:
        raise ValueError("Species must be unique and sample_per_species positive")
    out = output.resolve()
    if "raw" in out.parts or out.is_relative_to(root.resolve()):
        raise ValueError("Audit output must be outside raw/genome directories")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite prior results: {output}")
    source_hashes = {str(p): sha256(p) for p in (config_path, original, mapping)}
    manifest_ids, kingdom_counts, species_counts = {}, Counter(), Counter()
    accession_counts, fields_missing = Counter(), Counter()
    duplicate_ids = 0
    trailing_empty_rows = 0
    for row in rows(original, ("genome_id", "species_id", "fasta", "kingdom")):
        gid = row["genome_id"]
        trailing_empty_rows += row["_trailing_empty_fields"] > 0
        duplicate_ids += gid in manifest_ids
        manifest_ids[gid] = (row["species_id"], row["kingdom"])
        kingdom_counts[row["kingdom"]] += 1
        species_counts[(row["kingdom"], row["species_id"])] += 1
        if row["kingdom"] == "bacteria":
            accession_counts[gid.split("_", 1)[0]] += 1
        for field in ("genome_id", "species_id", "fasta", "kingdom"):
            value = row[field]
            fields_missing[field] += not bool(value)
    mapping_counts = Counter()
    selected = defaultdict(list)
    mapped_ids, target_paths, file_accessions = Counter(), Counter(), Counter()
    for row in rows(mapping, ("genome_id", "file_accession", "species_id", "kingdom", "target_relpath")):
        gid = row["genome_id"]
        mapping_counts["rows"] += 1
        mapped_ids[gid] += 1
        target_paths[row["target_relpath"]] += 1
        file_accessions[row["file_accession"]] += 1
        mapping_counts["genome_id_differs_from_file_accession"] += gid != row["file_accession"]
        mapping_counts["unknown_genome_id"] += gid not in manifest_ids
        mapping_counts["species_or_kingdom_mismatch"] += (
            gid in manifest_ids and manifest_ids[gid] != (row["species_id"], row["kingdom"])
        )
        if row["species_id"] in species and row["kingdom"] == "bacteria":
            selected[row["species_id"]].append(row)
    mapping_counts["duplicate_genome_id_excess_rows"] = sum(n - 1 for n in mapped_ids.values())
    mapping_counts["duplicate_target_path_excess_rows"] = sum(n - 1 for n in target_paths.values())
    mapping_counts["duplicate_file_accession_excess_rows"] = sum(n - 1 for n in file_accessions.values())
    inventory, samples = [], []
    for sp in species:
        candidates = selected[sp]
        directories = {str(Path(r["target_relpath"]).parent) for r in candidates}
        available = {}
        directory_errors = []
        for directory in sorted(directories):
            try:
                with os.scandir(safe_path(root, directory)) as entries:
                    for entry in entries:
                        if entry.name.endswith((".fna.gz", ".fna", ".fa.gz", ".fasta.gz")) and entry.is_file():
                            available[str(Path(directory) / entry.name)] = entry.stat().st_size
            except OSError as exc:
                directory_errors.append(f"{directory}: {exc}")
        targets = {r["target_relpath"] for r in candidates}
        present = [r for r in candidates if available.get(r["target_relpath"], 0) > 0]
        inventory.append({
            "species": sp,
            "manifest_records": species_counts[("bacteria", sp)],
            "mapped_records": len(candidates),
            "mapped_unique_genome_ids": len({r["genome_id"] for r in candidates}),
            "mapped_unique_paths": len(targets),
            "present_nonempty_records": len(present),
            "absent_or_empty_records": len(candidates) - len(present),
            "directory_fasta_files": len(available),
            "unmapped_directory_files": len(set(available) - targets),
            "genome_id_differs_from_file_accession": sum(r["genome_id"] != r["file_accession"] for r in candidates),
            "directory_errors": " | ".join(directory_errors),
        })
        # Select from mapped records BEFORE checking existence: no survivor substitution.
        chosen = sorted(candidates, key=lambda r: (sample_rank(config["seed"], r["genome_id"]), r["target_relpath"]))[:config["sample_per_species"]]
        for row in chosen:
            record = {k: row[k] for k in ("genome_id", "file_accession", "species_id", "target_relpath")}
            try:
                record.update(fasta_stats(safe_path(root, row["target_relpath"])))
                record["status"] = "ok"
                record["error"] = ""
            except (OSError, EOFError, ValueError, UnicodeError) as exc:
                record["status"] = "error"
                record["error"] = str(exc)
            samples.append(record)
    for path, digest in source_hashes.items():
        if sha256(path) != digest:
            raise ValueError(f"Input changed during audit: {path}")
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "species_inventory.tsv", inventory, list(inventory[0]))
    sample_columns = ["genome_id", "file_accession", "species_id", "target_relpath", "status", "error", "bytes", "file_sha256", "contig_multiset_sha256", "sequence_count", "total_bases", "n50_bases", "max_contig_bases", "n_fraction", "ambiguous_fraction", "gc_fraction_acgt"]
    write_tsv(output / "fasta_sample.tsv", samples, sample_columns)
    counts = [{"kingdom": k, "species": s, "records": n} for (k, s), n in sorted(species_counts.items())]
    write_tsv(output / "manifest_species_counts.tsv", counts, ["kingdom", "species", "records"])
    valid = [s for s in samples if s["status"] == "ok"]
    digests = Counter(s["contig_multiset_sha256"] for s in valid)
    summary = {
        "status": "completed_read_only_audit_not_training_readiness",
        "owner": "Codex current session",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "argv": sys.argv,
        "config": config,
        "input_sha256": source_hashes,
        "code_sha256": {str(p): sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
        "manifest": {"records": sum(kingdom_counts.values()), "duplicate_genome_id_excess_rows": duplicate_ids,
                     "rows_with_trailing_empty_columns": trailing_empty_rows,
                     "kingdom_counts": dict(kingdom_counts), "bacterial_species": sum(k == "bacteria" for k, s in species_counts),
                     "bacterial_accession_prefixes": dict(accession_counts), "blank_fields": dict(fields_missing)},
        "mapping": dict(mapping_counts),
        "sample": {"selection": "lowest SHA256(seed TAB genome_id), tie by target_relpath; selected before existence check",
                   "records": len(samples), "passed": len(valid), "failed": len(samples) - len(valid),
                   "exact_contig_multiset_duplicate_excess": sum(n - 1 for n in digests.values())},
        "limitations": ["File presence is not sequence validity; only sampled files were decoded and hashed.",
                        "No accession identity, biological contamination, annotation, license or clinical permissions verified.",
                        "File accession substitutions require upstream assembly metadata; matching numeric accession parts is insufficient.",
                        "No near-clone deduplication, lineage inference, labels or pretraining overlap assessment performed."]
    }
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    return summary
