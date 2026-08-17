#!/bin/bash
# PSO benchmark sweep for tab:pso-performance (paper appendix).
# Runs PSI-Card(3) / PSI-Sum(4) / PSI-Threshold(5) at n = 2^{16,18,20,22}
# over 10gbit / 1gbit / 100mbit / 10mbit through benchmark.sh, and appends
# one CSV line per run to log/benchmark/pso_sweep.csv:
#
#     pso,logn,bw,comm_MB,runtime_s
#
# Usage:  sudo ./pso_sweep.sh [bw ...]
#   e.g.  sudo ./pso_sweep.sh              # all four bandwidths (~1.5h)
#         sudo ./pso_sweep.sh 10gbit 1gbit # only these
#
# Individual run logs are archived as log/benchmark/pso_<idx>_<nn>_<bw>.*.log

set -u
cd "$(dirname "$(realpath "$0")")"

if [ "$(id -u)" -ne 0 ]; then echo "please run with sudo"; exit 1; fi

BWS=("$@")
[ ${#BWS[@]} -eq 0 ] && BWS=(10gbit 1gbit 100mbit 10mbit)

LOGDIR=log/benchmark
CSV=$LOGDIR/pso_sweep.csv
mkdir -p "$LOGDIR"
[ -f "$CSV" ] || echo "pso,logn,bw,comm_MB,runtime_s" > "$CSV"

declare -A NAME=([1]=ssPMT [2]=PSU [3]=PSI-Card [4]=PSI-Sum [5]=PSI-Threshold)

for bw in "${BWS[@]}"; do
  for idx in 1 2 3 4 5; do
    for nn in 16 18 20 22; do
      tag="${NAME[$idx]},$nn,$bw"
      # skip combos already recorded (lets the sweep resume after interruption)
      if grep "^$tag," "$CSV" | grep -qv FAIL; then echo "skip $tag (already done)"; continue; fi

      echo "=== ${NAME[$idx]}  n=2^$nn  $bw  ($(date +%H:%M:%S)) ==="
      ./benchmark.sh "$idx" -nn "$nn" -bw "$bw" -v > /dev/null 2>&1
      rc=$?
      rlog=$LOGDIR/last_recver.log
      cp "$rlog" "$LOGDIR/pso_${idx}_${nn}_${bw}.recver.log" 2>/dev/null
      cp "$LOGDIR/last_sender.log" "$LOGDIR/pso_${idx}_${nn}_${bw}.sender.log" 2>/dev/null

      clean=$(sed 's/\x1b\[[0-9;]*m//g' "$rlog")
      ms=$(echo "$clean" | grep -aoE "Passed +[0-9]+ms" | grep -oE "[0-9]+" | head -1)
      comm=$(echo "$clean" | grep -aoE "= [0-9.]+ MB" | head -1 | grep -oE "[0-9.]+")
      # ssPMT prints "recv A MB, sent B MB" instead
      [ -z "$comm" ] && comm=$(echo "$clean" | grep -aoE "recv [0-9.]+ MB, sent [0-9.]+ MB" \
             | head -1 | awk '{print $2 + $5}')
      if [ $rc -ne 0 ] || [ -z "$ms" ]; then
        echo "$tag,FAIL,FAIL" >> "$CSV"
        echo "    !! FAILED (rc=$rc) - see $LOGDIR/pso_${idx}_${nn}_${bw}.recver.log"
      else
        sec=$(awk -v m="$ms" 'BEGIN{printf "%.2f", m/1000}')
        echo "$tag,$comm,$sec" >> "$CSV"
        echo "    comm ${comm} MB, runtime ${sec} s"
      fi
    done
  done
done

chown "${SUDO_USER:-root}" "$CSV" "$LOGDIR"/pso_*.log 2>/dev/null
echo; echo "done. results in $CSV"
