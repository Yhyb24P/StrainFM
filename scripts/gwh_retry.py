#!/usr/bin/env python3
import concurrent.futures
import csv
import os
import subprocess
import time

ROOT = "/mnt/e/StrainFM"
GEN = ROOT + "/data/raw/genomes"
MAP = ROOT + "/data/raw/gwh_map_ecoli.tsv"
LOG = ROOT + "/logs/gwh_ecoli_retry.log"


def one(r):
    target = os.path.join(GEN, r["target_relpath"])
    if os.path.exists(target) and os.path.getsize(target) > 0:
        return "skip"
    os.makedirs(os.path.dirname(target), exist_ok=True)
    tmp = target + ".gwhretry"
    cmd = ["curl", "-sL", "--noproxy", "*", "--fail", "--retry", "5", "--retry-delay", "2",
           "--max-time", "300", "-o", tmp, r["gwh_url"]]
    x = subprocess.run(cmd, capture_output=True, text=True)
    if x.returncode != 0:
        try: os.remove(tmp)
        except FileNotFoundError: pass
        return f"FAIL {r['fasta_basename']} rc={x.returncode}"
    if subprocess.run(["gzip", "-t", tmp], capture_output=True).returncode != 0:
        os.remove(tmp)
        return f"GZIPFAIL {r['fasta_basename']}"
    os.replace(tmp, target)
    return "OK"


rows = list(csv.DictReader(open(MAP, newline=""), delimiter="\t"))
pending = [r for r in rows if not (os.path.exists(os.path.join(GEN, r["target_relpath"])) and os.path.getsize(os.path.join(GEN, r["target_relpath"])) > 0)]
print(f"pending={len(pending)}", flush=True)
ok = fail = 0
with open(LOG, "w") as log, concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
    for f in concurrent.futures.as_completed([ex.submit(one, r) for r in pending]):
        res = f.result()
        log.write(res + "\n"); log.flush()
        if res == "OK": ok += 1
        elif res.startswith("FAIL") or res.startswith("GZIPFAIL"): fail += 1
print(f"ok={ok} fail={fail}", flush=True)
