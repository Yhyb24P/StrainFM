#!/usr/bin/env bash
set -u

BASE=/mnt/e/StrainFM
RAW="$BASE/data/raw"
GEN="$RAW/genomes"
LOG="$BASE/logs/download_$(date +%Y%m%d_%H%M%S).log"
PROXY="http://127.0.0.1:28080"

exec >>"$LOG" 2>&1
echo "=== StrainFM download start $(date -Is) ==="
echo "log: $LOG"
echo "free before:"; df -h "$BASE" | tail -1

aria2_args=(
  -d "$GEN"
  -j 8
  -x 1
  -s 1
  -c
  --max-tries=5
  --retry-wait=3
  --timeout=120
  --connect-timeout=30
  --file-allocation=none
  --auto-file-renaming=false
  --allow-overwrite=false
  --console-log-level=warn
  --summary-interval=300
  --log-level=warn
  --all-proxy="$PROXY"
)

echo "=== E. coli start $(date -Is) ==="
aria2c -i "$RAW/ecoli_input.txt" "${aria2_args[@]}"
echo "=== E. coli end rc=$? $(date -Is) ==="
echo "ecoli files: $(find "$GEN/bacteria/Escherichia_coli" -type f -name '*.fna.gz' | wc -l)"
echo "ecoli size: $(du -sh "$GEN/bacteria/Escherichia_coli" | cut -f1)"

echo "=== remaining genomes start $(date -Is) ==="
aria2c -i "$RAW/others_input.txt" "${aria2_args[@]}"
echo "=== remaining genomes end rc=$? $(date -Is) ==="

echo "=== NCBI virus start $(date -Is) ==="
aria2c -i "$RAW/virus_input.txt" "${aria2_args[@]}"
echo "=== NCBI virus end rc=$? $(date -Is) ==="

echo "=== all files: $(find "$GEN" -type f -name '*.fna.gz' | wc -l) ==="
echo "=== total size: $(du -sh "$GEN" | cut -f1) ==="
echo "free after:"; df -h "$BASE" | tail -1
echo "=== StrainFM download end $(date -Is) ==="
