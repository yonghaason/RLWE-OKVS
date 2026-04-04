#!/usr/bin/env bash

set -Eeuo pipefail

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2
}

trap 'rc=$?; log "ERROR line=$LINENO cmd=${BASH_COMMAND} exit=$rc"; exit $rc' ERR

usage() {
  cat <<'EOF'
Usage:
  scripts/width_test_probe_fit.sh [options]

Options:
  --probe-nn <int>       target nn for width_test probe (default: 16)
  --source-nn <int>      source autosweep nn used for seed width/fit (default: 10)
  --ratio-start <float>  starting ratio (default: 2.0)
  --ratio-end <float>    ending ratio (default: 3.0)
  --ratio-step <float>   ratio step (default: 0.1)
  --th-cnt <int>         width_test thread count (default: 12)
  --test-index <int>     test index for OKVS_width_test (default: 0)
  --run-bin <path>       run binary path (default: ./build/run)
  --source-root <path>   autosweep root (default: log/okvs_width)
  --seed-summary <path>  optional TSV from prior probe; if given, use chosen_width as seed
  --out <path>           output TSV path
  --help                 show this help

Notes:
  - Probe runs use: -trials 100 -nosave
  - Width search target range is log2_failure_rate in [-9, -3]
  - If measured log2_failure_rate > -3, width is increased.
  - If measured log2_failure_rate < -9, width is decreased.
EOF
}

PROBE_NN=16
SOURCE_NN=10
RATIO_START=2.0
RATIO_END=3.0
RATIO_STEP=0.1
TH_CNT=12
TEST_INDEX=0
RUN_BIN=./build/run
SOURCE_ROOT=log/okvs_width
SEED_SUMMARY=""
OUT_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --probe-nn) PROBE_NN="$2"; shift 2 ;;
    --source-nn) SOURCE_NN="$2"; shift 2 ;;
    --ratio-start) RATIO_START="$2"; shift 2 ;;
    --ratio-end) RATIO_END="$2"; shift 2 ;;
    --ratio-step) RATIO_STEP="$2"; shift 2 ;;
    --th-cnt) TH_CNT="$2"; shift 2 ;;
    --test-index) TEST_INDEX="$2"; shift 2 ;;
    --run-bin) RUN_BIN="$2"; shift 2 ;;
    --source-root) SOURCE_ROOT="$2"; shift 2 ;;
    --seed-summary) SEED_SUMMARY="$2"; shift 2 ;;
    --out) OUT_FILE="$2"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

if [[ ! -x "$RUN_BIN" ]]; then
  echo "run binary not found: $RUN_BIN" >&2
  exit 1
fi

if [[ -z "$OUT_FILE" ]]; then
  OUT_FILE="log/okvs_nn${PROBE_NN}_probe_fit.tsv"
fi

mkdir -p "$(dirname "$OUT_FILE")"

find_source_summary() {
  local ratio_fixed="$1"
  find "$SOURCE_ROOT" -maxdepth 2 -type f \
    -path "*nn${SOURCE_NN}_ratio${ratio_fixed}_*/autosweep_summary.txt" \
    | sort | tail -n 1
}

seed_width_from_summary() {
  local ratio="$1"
  local seed_file="$2"
  python3 - "$ratio" "$seed_file" <<'PY'
import csv, sys
ratio = float(sys.argv[1])
path = sys.argv[2]
with open(path, newline='') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if abs(float(row['ratio']) - ratio) < 1e-9:
            print(row['chosen_width'])
            sys.exit(0)
sys.exit(1)
PY
}

parse_summary_field() {
  local file="$1"
  local key="$2"
  awk -F'=' -v k="$key" '$1 == k { print $2; exit }' "$file"
}

out_has_ratio() {
  local out_file="$1"
  local ratio="$2"
  [[ -f "$out_file" ]] || return 1
  awk -F'\t' -v r="$ratio" 'NR > 1 && $1 == r { found=1; exit } END { exit !found }' "$out_file"
}

parse_probe_metrics() {
  local file="$1"
  python3 - "$file" <<'PY'
from pathlib import Path
import re, sys
text = Path(sys.argv[1]).read_text(errors='ignore')
fr = re.search(r'failure_rate=([0-9.]+)', text)
lg = re.search(r'log2_failure_rate=([-a-zA-Z0-9.+]+)', text)
print((fr.group(1) if fr else ""), "\t", (lg.group(1) if lg else ""), sep="")
PY
}

ratio_values() {
  python3 - "$RATIO_START" "$RATIO_END" "$RATIO_STEP" <<'PY'
import sys
start = float(sys.argv[1])
end = float(sys.argv[2])
step = float(sys.argv[3])
vals = []
x = start
while x <= end + 1e-9:
    vals.append(f"{x:.2f}")
    x += step
print("\n".join(vals))
PY
}

if [[ ! -f "$OUT_FILE" ]]; then
  printf 'ratio\tseed_width\tchosen_width\tmeasured_failure_rate\tmeasured_log2_failure_rate\tfit_slope\tfit_intercept\textrapolated_w_at_log2_-40\tsource_summary\n' > "$OUT_FILE"
fi

while IFS= read -r ratio_fixed; do
  ratio_num="$(printf '%.2f' "$ratio_fixed")"
  if out_has_ratio "$OUT_FILE" "$ratio_num"; then
    log "skip ratio=${ratio_num}: already present in $OUT_FILE"
    continue
  fi

  log "start ratio=${ratio_num}"
  source_summary="$(find_source_summary "$ratio_fixed")"
  if [[ -z "$source_summary" ]]; then
    echo "Missing source summary for ratio=$ratio_fixed" >&2
    exit 1
  fi

  fit_slope="$(parse_summary_field "$source_summary" fit_slope)"

  if [[ -n "$SEED_SUMMARY" ]]; then
    seed_width="$(seed_width_from_summary "$ratio_fixed" "$SEED_SUMMARY")"
  else
    seed_width="$(parse_summary_field "$source_summary" full_stop_width)"
  fi

  log "ratio=${ratio_num} seed_width=${seed_width} slope=${fit_slope} source=${source_summary}"
  w="$seed_width"
  for _ in $(seq 1 64); do
    tmp_log="$(mktemp)"
    log "ratio=${ratio_num} run w=${w}"
    "$RUN_BIN" -u "$TEST_INDEX" -nn "$PROBE_NN" -ratio "$ratio_fixed" \
      -w0 "$w" -w1 "$((w + 1))" -th_cnt "$TH_CNT" -trials 100 -nosave > "$tmp_log"

    metrics="$(parse_probe_metrics "$tmp_log")"
    rm -f "$tmp_log"

    failure_rate="$(printf '%s' "$metrics" | cut -f1)"
    log2_failure_rate="$(printf '%s' "$metrics" | cut -f2)"
    log "ratio=${ratio_num} w=${w} failure_rate=${failure_rate:-NA} log2_failure_rate=${log2_failure_rate:-NA}"

    if [[ "$log2_failure_rate" == "-inf" || "$log2_failure_rate" == "inf" || "$log2_failure_rate" == "nan" || "$log2_failure_rate" == "NA" || -z "$log2_failure_rate" ]]; then
      w=$((w - 1))
      log "ratio=${ratio_num} adjust w->${w} because log2 metric is non-finite"
      continue
    fi

    if awk -v x="$log2_failure_rate" 'BEGIN { exit !(x > -3.0) }'; then
      w=$((w + 1))
      log "ratio=${ratio_num} adjust w->${w} because log2_failure_rate > -3"
      continue
    fi

    if awk -v x="$log2_failure_rate" 'BEGIN { exit !(x < -9.0) }'; then
      w=$((w - 1))
      log "ratio=${ratio_num} adjust w->${w} because log2_failure_rate < -9"
      continue
    fi

    fit_intercept="$(awk -v a="$fit_slope" -v x="$w" -v y="$log2_failure_rate" 'BEGIN { printf "%.16f", y - a * x }')"
    extrapolated="$(awk -v slope="$fit_slope" -v intercept="$fit_intercept" 'BEGIN { printf "%.6f", (-40.0 - intercept) / slope }')"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ratio_num" "$seed_width" "$w" "$failure_rate" "$log2_failure_rate" \
      "$fit_slope" "$fit_intercept" "$extrapolated" "$source_summary" >> "$OUT_FILE"
    log "ratio=${ratio_num} done chosen_width=${w} intercept=${fit_intercept} extrapolated=${extrapolated}"
    break
  done
done < <(ratio_values)

log "written: $OUT_FILE"
cat "$OUT_FILE"
