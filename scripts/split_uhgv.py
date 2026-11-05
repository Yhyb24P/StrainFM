#!/usr/bin/env python3
import csv
import gzip
import os
import re
import subprocess
import sys

manifest = "/mnt/c/Users/12087/Downloads/manifest_all.tsv"
input_zst = "/mnt/e/StrainFM/data/raw/uhgv/uhgv_full.fna.zst"
out_root = "/mnt/e/StrainFM/data/raw/genomes/virus"
log_path = "/mnt/e/StrainFM/logs/uhgv_split.log"

if len(sys.argv) >= 3:
    input_zst = sys.argv[1]
    out_root = sys.argv[2]

mapping = {}
with open(manifest, newline="", encoding="utf-8", errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        gid = (row.get("genome_id") or "").strip()
        if gid.startswith("UHGV-"):
            mapping[gid] = (row.get("species_id") or "").strip()

os.makedirs(out_root, exist_ok=True)
count = 0
written = set()
missing = []
current = None
out_handle = None
header_re = re.compile(r"UHGV-\d+")

def close_current():
    global out_handle
    if out_handle is not None:
        out_handle.close()
        out_handle = None

with open(log_path, "w", encoding="utf-8") as log:
    log.write(f"input={input_zst}\nmanifest_uhgv={len(mapping)}\n")
    proc = subprocess.Popen(["zstd", "--long=29", "--memory=1024MB", "-dc", input_zst], stdout=subprocess.PIPE, text=True, errors="replace")
    for line in proc.stdout:
        if line.startswith(">"):
            close_current()
            m = header_re.search(line)
            if not m:
                current = None
                log.write(f"WARN unparsed header: {line[:120]!r}\n")
                continue
            current = m.group(0)
            species = mapping.get(current)
            if species is None:
                missing.append(current)
                log.write(f"WARN not in manifest: {current}\n")
                current = None
                continue
            target_dir = os.path.join(out_root, species)
            os.makedirs(target_dir, exist_ok=True)
            target = os.path.join(target_dir, current + ".fna.gz")
            out_handle = gzip.open(target, "wt", encoding="utf-8")
            out_handle.write(line)
            count += 1
            written.add(current)
        elif current is not None and out_handle is not None:
            out_handle.write(line)
    close_current()
    rc = proc.wait()
    log.write(f"zstd_rc={rc}\nrecords={count}\nunique={len(written)}\nmanifest_missing={len(set(mapping)-written)}\nheader_missing={len(missing)}\n")
    for gid in sorted(set(mapping) - written)[:50]:
        log.write(f"MISSING {gid}\n")
print(f"records={count} unique={len(written)} missing={len(set(mapping)-written)}")
