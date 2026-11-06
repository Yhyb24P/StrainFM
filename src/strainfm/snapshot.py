"""Resumable full-species input snapshot. The source files are never modified."""

import json
import os
import platform
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from .audit import fasta_stats, rows, safe_path, sha256, write_tsv


def inspect_file(item):
    record, path = item
    result = dict(record)
    try:
        before = path.stat()
        result.update(observed_size=before.st_size, observed_mtime_ns=before.st_mtime_ns)
        if (before.st_size, before.st_mtime_ns) != (record["initial_size"], record["initial_mtime_ns"]):
            raise ValueError("File changed since inventory capture")
        if not before.st_size:
            raise ValueError("Empty input file")
        result.update(fasta_stats(path))
        result.update(status="valid", error="")
    except (OSError, EOFError, ValueError, UnicodeError) as exc:
        result.update(status="missing" if isinstance(exc, FileNotFoundError) else "invalid", error=str(exc))
    return result


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
    temporary.replace(path)


def prepare_output(output, contract):
    output.mkdir(parents=True, exist_ok=True)
    path = output / "contract.json"
    if path.exists():
        if json.loads(path.read_text()) != contract:
            raise ValueError("Resume contract changed; use a new output directory")
    elif list(output.iterdir()):
        raise ValueError("Refusing unrelated nonempty output directory")
    else:
        atomic_json(path, contract)


def run_snapshot(config_path, output, workers=8):
    config_path, output = Path(config_path), Path(output)
    config = json.loads(config_path.read_text())
    mapping, root = Path(config["download_manifest"]), Path(config["genome_root"])
    if "raw" in output.resolve().parts or output.resolve().is_relative_to(root.resolve()):
        raise ValueError("Output must be outside raw/genome directories")
    if not 1 <= workers <= 24:
        raise ValueError("workers must be between 1 and 24")
    code = {str(p): sha256(p) for p in (Path(__file__), Path(__file__).with_name("audit.py"))}
    contract = {"schema": 1, "config": config, "config_sha256": sha256(config_path), "mapping_sha256": sha256(mapping), "code_sha256": code}
    prepare_output(output, contract)
    contract_path = output / "contract.json"
    inventory_path = output / "inventory.json"
    if not inventory_path.exists():
        selected = [{k: r[k] for k in ("genome_id", "file_accession", "species_id", "target_relpath")} for r in rows(mapping, ("genome_id", "file_accession", "species_id", "kingdom", "target_relpath")) if r["species_id"] == config["species"] and r["kingdom"] == "bacteria"]
        if not selected:
            raise ValueError("No matching bacterial assemblies")
        if len({r["genome_id"] for r in selected}) != len(selected) or len({r["target_relpath"] for r in selected}) != len(selected):
            raise ValueError("Non-unique genome IDs or file paths")
        by_directory = {}
        for relative in {str(Path(r["target_relpath"]).parent) for r in selected}:
            directory = safe_path(root, relative)
            try:
                with os.scandir(directory) as entries:
                    by_directory[relative] = {entry.name: entry.stat() for entry in entries if entry.is_file()}
            except FileNotFoundError:
                by_directory[relative] = {}
        for r in selected:
            rel = Path(r["target_relpath"])
            if rel.is_absolute() or ".." in rel.parts:
                raise ValueError("Unsafe input path")
            stat = by_directory[str(rel.parent)].get(rel.name)
            r.update(initial_size=stat.st_size if stat else -1, initial_mtime_ns=stat.st_mtime_ns if stat else -1)
        atomic_json(inventory_path, sorted(selected, key=lambda r: r["genome_id"]))
    inventory = json.loads(inventory_path.read_text())
    journal_path = output / "file_checks.jsonl"
    completed = {}
    if journal_path.exists():
        # A crash may leave only the final journal line incomplete.
        with journal_path.open("rb+") as handle:
            while True:
                start = handle.tell()
                line = handle.readline()
                if not line:
                    break
                if not line.endswith(b"\n"):
                    handle.truncate(start)
                    break
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    if handle.read():
                        raise ValueError("Corruption inside snapshot journal")
                    handle.truncate(start)
                    break
                if row["genome_id"] in completed:
                    raise ValueError("Duplicate snapshot journal record")
                completed[row["genome_id"]] = row
    pending = [(r, root / r["target_relpath"]) for r in inventory if r["genome_id"] not in completed]
    started = datetime.now(timezone.utc).isoformat()
    atomic_json(output / "progress.json", {"status": "running", "pid": os.getpid(), "started_utc": started, "completed": len(completed), "total": len(inventory)})
    with journal_path.open("a") as journal, ProcessPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(inspect_file, pending, chunksize=8):
            journal.write(json.dumps(result, ensure_ascii=False) + "\n")
            completed[result["genome_id"]] = result
            if len(completed) % 250 == 0:
                journal.flush()
                atomic_json(output / "progress.json", {"status": "running", "pid": os.getpid(), "completed": len(completed), "total": len(inventory), "updated_utc": datetime.now(timezone.utc).isoformat()})
                print(f"snapshot {len(completed)}/{len(inventory)}", flush=True)
    if sha256(mapping) != contract["mapping_sha256"]:
        raise ValueError("Input mapping changed during snapshot")
    records = [completed[r["genome_id"]] for r in inventory]
    columns = list(dict.fromkeys(k for r in records for k in r))
    write_tsv(output / "assemblies.tsv", records, columns)
    valid = [r for r in records if r["status"] == "valid"]
    digests = Counter(r["contig_multiset_sha256"] for r in valid)
    result = {"status": "file_snapshot_complete", "species": config["species"], "finished_utc": datetime.now(timezone.utc).isoformat(), "records": len(records), "file_status": dict(Counter(r["status"] for r in records)), "exact_record_multiset_groups": len(digests), "exact_record_multiset_duplicate_excess": sum(n - 1 for n in digests.values()), "compressed_bytes": sum(r["bytes"] for r in valid), "assemblies_sha256": sha256(output / "assemblies.tsv"), "contract_sha256": sha256(contract_path), "inventory_sha256": sha256(inventory_path), "python": platform.python_version(), "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(), "limitations": ["Exact FASTA record multiset digest does not normalize reverse complements, circular starts or assembly breakpoints.", "File validation is not species identity, biological quality or epidemiological independence confirmation.", "Source files remain in place; later use must verify content hashes against this snapshot."]}
    atomic_json(output / "summary.json", result)
    atomic_json(output / "progress.json", {"status": "completed", "completed": len(records), "total": len(records)})
    return result
