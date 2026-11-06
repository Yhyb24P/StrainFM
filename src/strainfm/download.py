"""Download and validate genome FASTA files before replacing a destination."""

import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .audit import fasta_stats, rows, safe_path


def valid_fasta(path):
    try:
        fasta_stats(path)
        return True
    except (OSError, EOFError, ValueError, UnicodeError):
        return False


def download_fasta(root, relative, url):
    target = safe_path(Path(root), relative)
    if valid_fasta(target):
        return "skip"
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target.parent, prefix=target.name + ".part-", suffix=target.suffix, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        subprocess.run(["curl", "-sSL", "--fail", "--retry", "5", "--retry-delay", "3", "--max-time", "300", "-o", str(temporary), url], check=True, capture_output=True)
        fasta_stats(temporary)
        temporary.replace(target)
        return "OK"
    except (OSError, EOFError, ValueError, UnicodeError, subprocess.CalledProcessError) as exc:
        return f"FAIL {relative}: {exc}"
    finally:
        temporary.unlink(missing_ok=True)


def run_download(manifest, root):
    records = list(rows(manifest, ("target_relpath",)))
    for row in records:
        safe_path(Path(root), row["target_relpath"])
        row["url"] = row.get("gwh_url") or row.get("ncbi_url")
        if not row["url"]:
            raise ValueError(f"Missing download URL: {row['target_relpath']}")
    with ThreadPoolExecutor(max_workers=8) as pool:
        counts = Counter()
        for result in pool.map(lambda row: download_fasta(root, row["target_relpath"], row["url"]), records):
            counts["failed" if result.startswith("FAIL") else result] += 1
            if result.startswith("FAIL"):
                print(result, flush=True)
    return {"records": len(records), "downloaded": counts["OK"], "skipped_valid": counts["skip"], "failed": counts["failed"]}
