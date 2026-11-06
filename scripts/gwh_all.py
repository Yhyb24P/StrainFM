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
from collections import defaultdict
from urllib.parse import quote

from strainfm.download import download_fasta, valid_fasta

ROOT = "/mnt/e/StrainFM"
GEN = ROOT + "/data/raw/genomes"
API = "https://ngdc.cncb.ac.cn/gwh/api/esdt/gwh_assemblies/query"
BASE = "https://download.cncb.ac.cn/assembly/ncbi/"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
LOG = ROOT + "/logs/gwh_all_" + time.strftime("%Y%m%d_%H%M%S") + ".log"
MAP = ROOT + "/data/raw/gwh_map_all.tsv"
MISS = ROOT + "/data/raw/gwh_misses_all.tsv"
NAMES = ["scientificName", "commonNames", "group", "source", "accession",
         "refseqAccession", "hasReannotation", "genomeRepresentation", "assemblyLevel",
         "genomeSize", "chrCount", "gcContent", "releaseDate", "sourceReleaseDate",
         "dna", "gff", "rna", "cds", "protein", "refseqDna", "refseqGff",
         "refseqRna", "refseqCds", "refseqProtein"]

lock = threading.Lock()
logf = open(LOG, "w")
stats = defaultdict(int)


def log(msg):
    with lock:
        line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {msg}"
        logf.write(line + "\n"); logf.flush()
        print(line, flush=True)


def api_page(species, start):
    cols = [{"data": "", "name": "", "searchable": False, "orderable": False,
             "search": {"value": "", "regex": False}}]
    for n in NAMES:
        cols.append({"data": n, "name": n, "searchable": True, "orderable": True,
                     "search": {"value": species if n == "scientificName" else "", "regex": False}})
    payload = {"draw": 1, "columns": cols, "order": [{"column": 1, "dir": "asc"}],
               "start": start, "length": 1000, "search": {"value": "", "regex": False}}
    for attempt in range(3):
        try:
            r = subprocess.run(["curl", "-sS", "--max-time", "60", "-A", UA, "-X", "POST", API,
                                "-H", "Content-Type: application/json", "--data", json.dumps(payload)],
                               capture_output=True, text=True)
            return json.loads(r.stdout)
        except Exception as e:
            time.sleep(1 + attempt)
    return None

def main():
    log("indexing existing files")
    existing = set()
    for dirpath, _, files in os.walk(GEN):
        for f in files:
            if f.endswith(".fna.gz") and valid_fasta(os.path.join(dirpath, f)):
                existing.add(os.path.relpath(os.path.join(dirpath, f), GEN))
    log(f"existing={len(existing)}")
    # Build remaining rows from both manifests
    remaining = []
    for src in [ROOT + "/data/raw/download_manifest_all.tsv",
                ROOT + "/data/raw/virus_download_manifest.tsv"]:
        if not os.path.exists(src):
            continue
        for r in csv.DictReader(open(src, newline=""), delimiter="\t"):
            if r["target_relpath"] not in existing:
                remaining.append(r)
    log(f"remaining={len(remaining)}")
    by_species = defaultdict(list)
    for r in remaining:
        by_species[r["species_id"]].append(r)
    species_sorted = sorted(by_species.items(), key=lambda kv: len(kv[1]), reverse=True)
    log(f"species={len(species_sorted)}")
    map_f = open(MAP, "w", newline="")
    miss_f = open(MISS, "w", newline="")
    mw = csv.DictWriter(map_f, delimiter="\t", fieldnames=["file_accession", "fasta_basename", "target_relpath", "species_id", "ftpDir", "dna", "gwh_url", "ncbi_url"])
    missw = csv.DictWriter(miss_f, delimiter="\t", fieldnames=["file_accession", "fasta_basename", "target_relpath", "species_id", "ncbi_url"])
    mw.writeheader(); missw.writeheader()
    for idx, (species, rows) in enumerate(species_sorted, 1):
        accs = {r["file_accession"]: r for r in rows}
        search = species.replace("_", " ")
        hits = {}
        start = 0
        total = None
        while True:
            d = api_page(search, start)
            if d is None:
                log(f"API_FAIL species={species} start={start}")
                break
            total = d.get("recordsFiltered", 0)
            data = d.get("data", [])
            for row in data:
                acc = row.get("accession")
                if acc in accs and row.get("source") == "NCBI" and row.get("ftpDir") and row.get("dna"):
                    r = accs[acc]
                    url = BASE + quote(row["ftpDir"] + "/" + row["dna"], safe="/._-~")
                    hits[acc] = {"file_accession": acc, "fasta_basename": r["fasta_basename"],
                                 "target_relpath": r["target_relpath"], "species_id": species,
                                 "ftpDir": row["ftpDir"], "dna": row["dna"], "gwh_url": url,
                                 "ncbi_url": r.get("ncbi_url", "")}
            if not data or start + len(data) >= total:
                break
            start += len(data)
            time.sleep(0.2)
        # download hits for this species
        ok = skip = fail = 0
        if hits:
            with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
                for f in concurrent.futures.as_completed([ex.submit(download_fasta, GEN, rec["target_relpath"], rec["gwh_url"]) for rec in hits.values()]):
                    res = f.result()
                    if res == "OK": ok += 1
                    elif res == "skip": skip += 1
                    else: fail += 1
                    if res.startswith("FAIL"):
                        log(res)
            for rec in hits.values():
                mw.writerow(rec); map_f.flush()
                stats["mapped"] += 1
            stats["downloaded"] += ok; stats["skipped"] += skip; stats["failed"] += fail
        # misses are rows with no GWH mapping
        for acc, r in accs.items():
            if acc not in hits:
                missw.writerow({k: r.get(k, "") for k in missw.fieldnames})
        miss_f.flush()
        log(f"[{idx}/{len(species_sorted)}] {species} rows={len(rows)} hit={len(hits)} ok={ok} skip={skip} fail={fail} gwh_total={total}")
    log(f"done {dict(stats)}")
    map_f.close(); miss_f.close(); logf.close()


if __name__ == "__main__":
    main()
