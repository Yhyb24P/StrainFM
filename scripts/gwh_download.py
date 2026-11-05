#!/usr/bin/env python3
# EXPERIMENTAL: uses the GWH browse-page DataTables endpoint, which is undocumented.
# Not a public/stable API and not part of the main pipeline; see DECISIONS.md D0-12.
import concurrent.futures
import csv
import json
import os
import subprocess
import sys
import threading
import time
from urllib.parse import quote

ROOT = "/mnt/e/StrainFM"
GENOMES = ROOT + "/data/raw/genomes"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
API = "https://ngdc.cncb.ac.cn/gwh/api/esdt/gwh_assemblies/query"
BASE = "https://download.cncb.ac.cn/assembly/ncbi/"
NAMES = [
    "scientificName", "commonNames", "group", "source", "accession",
    "refseqAccession", "hasReannotation", "genomeRepresentation", "assemblyLevel",
    "genomeSize", "chrCount", "gcContent", "releaseDate", "sourceReleaseDate",
    "dna", "gff", "rna", "cds", "protein", "refseqDna", "refseqGff",
    "refseqRna", "refseqCds", "refseqProtein",
]

lock = threading.Lock()
stats = {"hit": 0, "miss": 0, "downloaded": 0, "download_failed": 0, "skipped": 0}
log_fh = None


def log(msg):
    with lock:
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        print(line)
        log_fh.write(line + "\n")
        log_fh.flush()


def query_gwh(acc):
    cols = [{"data": "", "name": "", "searchable": False, "orderable": False,
             "search": {"value": "", "regex": False}}]
    for n in NAMES:
        cols.append({"data": n, "name": n, "searchable": True, "orderable": True,
                     "search": {"value": acc if n == "accession" else "", "regex": False}})
    payload = {"draw": 1, "columns": cols, "order": [{"column": 1, "dir": "asc"}],
               "start": 0, "length": 1, "search": {"value": "", "regex": False}}
    for attempt in range(3):
        try:
            r = subprocess.run(
                ["curl", "-sS", "--max-time", "25", "-A", UA, "-X", "POST", API,
                 "-H", "Content-Type: application/json", "--data", json.dumps(payload)],
                capture_output=True, text=True)
            d = json.loads(r.stdout)
            rows = d.get("data", [])
            if rows and rows[0].get("accession") == acc:
                row = rows[0]
                if row.get("ftpDir") and row.get("dna") and row.get("source") == "NCBI":
                    return row
            return None
        except Exception:
            time.sleep(1 + attempt)
    return None


def download_one(row, target_rel, basename, gwh_url):
    target = os.path.join(GENOMES, target_rel)
    os.makedirs(os.path.dirname(target), exist_ok=True)
    if os.path.exists(target) and os.path.getsize(target) > 0:
        with lock:
            stats["skipped"] += 1
        return "skipped"
    tmp = target + ".gwhpart"
    cmd = ["curl", "-sL", "--noproxy", "*", "--fail", "--retry", "3", "--retry-delay", "2",
           "--max-time", "300", "-o", tmp, gwh_url]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        try: os.remove(tmp)
        except FileNotFoundError: pass
        with lock:
            stats["download_failed"] += 1
        log(f"DOWNLOAD_FAIL {basename} rc={r.returncode} {r.stderr.strip()[:120]}")
        return "failed"
    if subprocess.run(["gzip", "-t", tmp], capture_output=True).returncode != 0:
        os.remove(tmp)
        with lock:
            stats["download_failed"] += 1
        log(f"GZIP_FAIL {basename}")
        return "failed"
    os.replace(tmp, target)
    with lock:
        stats["downloaded"] += 1
        n = stats["downloaded"]
    if n % 100 == 0:
        log(f"downloaded={n} hit={stats['hit']} miss={stats['miss']} failed={stats['download_failed']}")
    return "ok"


def main():
    global log_fh
    if len(sys.argv) < 3:
        print("usage: gwh_download.py INPUT.tsv PREFIX")
        sys.exit(2)
    inp, prefix = sys.argv[1], sys.argv[2]
    os.makedirs(ROOT + "/logs", exist_ok=True)
    log_fh = open(f"{ROOT}/logs/gwh_{prefix}_{time.strftime('%Y%m%d_%H%M%S')}.log", "w")
    rows = list(csv.DictReader(open(inp, newline=""), delimiter="\t"))
    log(f"start input={inp} rows={len(rows)}")
    mapping = []
    misses = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(query_gwh, r["file_accession"]): r for r in rows}
        for f in concurrent.futures.as_completed(futs):
            r = futs[f]
            row = f.result()
            if row is None:
                misses.append(r)
                with lock: stats["miss"] += 1
                continue
            gwh_url = BASE + quote(row["ftpDir"] + "/" + row["dna"], safe="/._-~")
            mapping.append({
                "file_accession": r["file_accession"],
                "fasta_basename": r["fasta_basename"],
                "target_relpath": r["target_relpath"],
                "ftpDir": row["ftpDir"],
                "dna": row["dna"],
                "gwh_url": gwh_url,
                "ncbi_url": r.get("ncbi_url", ""),
            })
            with lock: stats["hit"] += 1
    log(f"mapping done hit={stats['hit']} miss={stats['miss']}")
    map_path = f"{ROOT}/data/raw/gwh_map_{prefix}.tsv"
    with open(map_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(mapping[0].keys()) if mapping else ["file_accession"])
        w.writeheader(); w.writerows(mapping)
    miss_path = f"{ROOT}/data/raw/gwh_misses_{prefix}.tsv"
    with open(miss_path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=["file_accession", "fasta_basename", "target_relpath", "ncbi_url"])
        w.writeheader(); w.writerows({k: r.get(k, "") for k in w.fieldnames} for r in misses)
    log(f"mapping saved {map_path} misses {miss_path}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
        futs = [ex.submit(download_one, r["file_accession"], r["target_relpath"], r["fasta_basename"], r["gwh_url"]) for r in mapping]
        for f in concurrent.futures.as_completed(futs):
            f.result()
    log(f"done {stats}")
    log_fh.close()


if __name__ == "__main__":
    main()
