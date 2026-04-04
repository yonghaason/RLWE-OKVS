#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/width_test_autosweep.sh --ratio <bandExpansion> [options]

Required:
  --ratio <float>          bandExpansion ratio for width_test

Options:
  --run-bin <path>         executable path (default: ./build/run)
  --nn <int>               log2(n) input for width_test (default: 10)
  --logp <int>             plaintext modulus logp (default: 60)
  --start-w <int>          starting width for coarse search (default: 1)
  --max-w <int>            maximum width to try (default: 256)
  --th-cnt <int>           thread count (default: 12)

  --coarse-trials <int>    trials for coarse search (default: 100)
  --full-trials <int>      trials for full sweep (default: 100000)
  --coarse-min-fr <float>  lower bound of acceptable coarse failure rate (default: 0.10)
  --coarse-max-fr <float>  upper bound of acceptable coarse failure rate (default: 0.30)

  --target-log2-fr <float> stop full sweep once log2(failure_rate) < this value (default: -9)
  --target-fr <float>      explicit stop threshold on failure_rate (overrides target-log2-fr)

  --out-root <path>        root directory for autosweep outputs (default: tuning/okvs_width)
  --tag <string>           optional tag added to run directory name
  --extra-arg <token>      extra single token passed through to ./build/run (repeatable)

Persistent outputs per run:
  - autosweep_summary.txt
  - full_search.tsv

width_test itself will still write only its TSV under tuning/okvs_width/.
EOF
}

RUN_BIN="./build/run"
NN=10
RATIO=""
LOGP=60
START_W=1
MAX_W=256
TH_CNT=12

COARSE_TRIALS=100
FULL_TRIALS=100000
COARSE_MIN_FR=0.01
COARSE_MAX_FR=0.90

TARGET_LOG2_FR=-9
TARGET_FR=""

OUT_ROOT="tuning/okvs_width"
TAG=""
declare -a EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-bin) RUN_BIN="$2"; shift 2 ;;
    --nn) NN="$2"; shift 2 ;;
    --ratio) RATIO="$2"; shift 2 ;;
    --logp) LOGP="$2"; shift 2 ;;
    --start-w) START_W="$2"; shift 2 ;;
    --max-w) MAX_W="$2"; shift 2 ;;
    --th-cnt) TH_CNT="$2"; shift 2 ;;
    --coarse-trials) COARSE_TRIALS="$2"; shift 2 ;;
    --full-trials) FULL_TRIALS="$2"; shift 2 ;;
    --coarse-min-fr) COARSE_MIN_FR="$2"; shift 2 ;;
    --coarse-max-fr) COARSE_MAX_FR="$2"; shift 2 ;;
    --target-log2-fr) TARGET_LOG2_FR="$2"; shift 2 ;;
    --target-fr) TARGET_FR="$2"; shift 2 ;;
    --out-root) OUT_ROOT="$2"; shift 2 ;;
    --tag) TAG="$2"; shift 2 ;;
    --extra-arg) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "[ERROR] unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ -z "$RATIO" ]]; then
  echo "[ERROR] --ratio is required" >&2
  usage >&2
  exit 1
fi

if [[ ! -x "$RUN_BIN" ]]; then
  echo "[ERROR] executable not found or not executable: $RUN_BIN" >&2
  exit 1
fi

if (( START_W < 0 || MAX_W < START_W )); then
  echo "[ERROR] invalid width range: start=$START_W max=$MAX_W" >&2
  exit 1
fi

if [[ -z "$TARGET_FR" ]]; then
  TARGET_FR="$(awk -v x="$TARGET_LOG2_FR" 'BEGIN { printf "%.17g\n", exp(log(2)*x) }')"
fi

float_ge() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a >= b) }'; }
float_le() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a <= b) }'; }
float_lt() { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a < b) }'; }

log2_or_na() {
  awk -v x="$1" 'BEGIN {
    if (x > 0) printf "%.17g\n", log(x)/log(2);
    else print "NA";
  }'
}

format_rate_pow2() {
  local fr="$1"
  local lg2
  lg2="$(log2_or_na "$fr")"
  if [[ "$lg2" == "NA" ]]; then
    printf '0'
  else
    awk -v x="$lg2" 'BEGIN { printf "2^%.4f", x }'
  fi
}

parse_failure_rate() {
  local file="$1"
  local fr

  fr="$(awk -F',' '
    /^[[:space:]]*->[[:space:]]*failures=/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /^[[:space:]]*failure_rate[[:space:]]*=/) {
          sub(/^[[:space:]]*failure_rate[[:space:]]*=[[:space:]]*/, "", $i)
          val = $i
        }
      }
    }
    END {
      if (val != "") print val
      else exit 1
    }
  ' "$file" 2>/dev/null || true)"
  if [[ -n "$fr" ]]; then
    printf '%s\n' "$fr"
    return 0
  fi

  fr="$(awk -F'|' '
    /Failure[[:space:]_][Rr]ate:/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /Failure[[:space:]_][Rr]ate:/) {
          sub(/.*Failure[[:space:]_][Rr]ate:[[:space:]]*/, "", $i)
          sub(/[[:space:]]*=.*/, "", $i)
          gsub(/[[:space:]]+$/, "", $i)
          val = $i
        }
      }
    }
    END {
      if (val != "") print val
      else exit 1
    }
  ' "$file" 2>/dev/null || true)"
  if [[ -n "$fr" ]]; then
    printf '%s\n' "$fr"
    return 0
  fi

  return 1
}

parse_failures_pair() {
  local file="$1"
  local pair

  pair="$(awk '
    /^[[:space:]]*->[[:space:]]*failures=/ {
      if (match($0, /failures=[0-9]+\/[0-9]+/)) {
        s = substr($0, RSTART, RLENGTH)
        sub(/^failures=/, "", s)
        gsub(/\//, "\t", s)
        val = s
      }
    }
    END {
      if (val != "") print val
      else exit 1
    }
  ' "$file" 2>/dev/null || true)"
  if [[ -n "$pair" ]]; then
    printf '%s\n' "$pair"
    return 0
  fi

  pair="$(awk '
    /Failures:[[:space:]]*[0-9]+[[:space:]]*\/[[:space:]]*[0-9]+/ {
      if (match($0, /Failures:[[:space:]]*[0-9]+[[:space:]]*\/[[:space:]]*[0-9]+/)) {
        s = substr($0, RSTART, RLENGTH)
        sub(/^Failures:[[:space:]]*/, "", s)
        gsub(/[[:space:]]*\/[[:space:]]*/, "\t", s)
        val = s
      }
    }
    END {
      if (val != "") print val
      else exit 1
    }
  ' "$file" 2>/dev/null || true)"
  if [[ -n "$pair" ]]; then
    printf '%s\n' "$pair"
    return 0
  fi

  return 1
}

parse_avg_ms() {
  local file="$1"
  local ms

  ms="$(awk -F',' '
    /^[[:space:]]*->[[:space:]]*failures=/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /^[[:space:]]*avg_ms[[:space:]]*=/) {
          sub(/^[[:space:]]*avg_ms[[:space:]]*=[[:space:]]*/, "", $i)
          val = $i
        }
      }
    }
    END {
      if (val != "") print val
      else exit 1
    }
  ' "$file" 2>/dev/null || true)"
  if [[ -n "$ms" ]]; then
    printf '%s\n' "$ms"
    return 0
  fi

  ms="$(awk -F'|' '
    /Avg[[:space:]]+Time[[:space:]]+per[[:space:]]+Trial:/ {
      for (i = 1; i <= NF; ++i) {
        if ($i ~ /Avg[[:space:]]+Time[[:space:]]+per[[:space:]]+Trial:/) {
          sub(/.*Avg[[:space:]]+Time[[:space:]]+per[[:space:]]+Trial:[[:space:]]*/, "", $i)
          sub(/[[:space:]]*ms.*/, "", $i)
          gsub(/[[:space:]]+$/, "", $i)
          val = $i
        }
      }
    }
    END {
      if (val != "") print val
      else exit 1
    }
  ' "$file" 2>/dev/null || true)"
  if [[ -n "$ms" ]]; then
    printf '%s\n' "$ms"
    return 0
  fi

  return 1
}

ratio_fixed="$(awk -v r="$RATIO" 'BEGIN { printf "%.2f", r }')"
timestamp="$(date +"%y%m%d_%H%M%S")"
run_name="nn${NN}_ratio${ratio_fixed}_${timestamp}"
if [[ -n "$TAG" ]]; then
  run_name="${run_name}_${TAG}"
fi

RUN_DIR="${OUT_ROOT}/${run_name}"
mkdir -p "$RUN_DIR"

AUTO_FULL_TSV="${RUN_DIR}/full_search.tsv"
AUTO_SUMMARY="${RUN_DIR}/autosweep_summary.txt"

cat > "$AUTO_FULL_TSV" <<'EOF'
width	failure_rate	log2_failure_rate	failures	trials	avg_ms	status
EOF

{
  echo "run_bin=$RUN_BIN"
  echo "nn=$NN"
  echo "ratio=$RATIO"
  echo "ratio_fixed=$ratio_fixed"
  echo "logp=$LOGP"
  echo "start_w=$START_W"
  echo "max_w=$MAX_W"
  echo "th_cnt=$TH_CNT"
  echo "coarse_trials=$COARSE_TRIALS"
  echo "full_trials=$FULL_TRIALS"
  echo "coarse_min_fr=$COARSE_MIN_FR"
  echo "coarse_max_fr=$COARSE_MAX_FR"
  echo "target_log2_fr=$TARGET_LOG2_FR"
  echo "target_fr=$TARGET_FR"
  echo "run_dir=$RUN_DIR"
} > "$AUTO_SUMMARY"

TMP_LOG="$(mktemp)"
cleanup() {
  rm -f "$TMP_LOG"
}
trap cleanup EXIT

run_one() {
  local stage="$1"
  local w="$2"
  shift 2

  "$@" > "$TMP_LOG" 2>&1

  local fr fail_pair avg_ms log2_fr failures trials
  fr="$(parse_failure_rate "$TMP_LOG")" || {
    echo "[ERROR] failed to parse failure_rate for stage=$stage w=$w" >&2
    cat "$TMP_LOG" >&2
    exit 3
  }
  fail_pair="$(parse_failures_pair "$TMP_LOG" || true)"
  failures="$(printf '%s' "$fail_pair" | cut -f1)"
  trials="$(printf '%s' "$fail_pair" | cut -f2)"
  avg_ms="$(parse_avg_ms "$TMP_LOG" || printf 'NA\n')"
  log2_fr="$(log2_or_na "$fr")"

  printf '%s\t%s\t%s\t%s\t%s\n' "$fr" "$log2_fr" "${failures:-NA}" "${trials:-NA}" "$avg_ms"
}

echo "[autosweep] outputs: $RUN_DIR"

FOUND_W=""
FOUND_FR=""
for (( w=START_W; w<=MAX_W; ++w )); do
  cmd=(
    "$RUN_BIN"
    -u 0
    -nn "$NN"
    -ratio "$RATIO"
    -logp "$LOGP"
    -w0 "$w"
    -w1 "$((w+1))"
    -th_cnt "$TH_CNT"
    -trials "$COARSE_TRIALS"
  )
  if ((${#EXTRA_ARGS[@]})); then
    cmd+=("${EXTRA_ARGS[@]}")
  fi

  result="$(run_one coarse "$w" "${cmd[@]}")"
  IFS=$'\t' read -r fr log2_fr failures trials avg_ms <<< "$result"

  if float_ge "$fr" "$COARSE_MIN_FR" && float_le "$fr" "$COARSE_MAX_FR"; then
    FOUND_W="$w"
    FOUND_FR="$fr"
    echo "[coarse] selected w=$w, failure_rate=$(printf '%.6f' "$fr"), ~$(format_rate_pow2 "$fr")"
    break
  fi
done

if [[ -z "$FOUND_W" ]]; then
  {
    echo "coarse_result=not_found"
    echo "reason=no width in coarse range produced actual failure_rate within [${COARSE_MIN_FR}, ${COARSE_MAX_FR}]"
  } >> "$AUTO_SUMMARY"

  echo "[ERROR] coarse search failed: no width met the actual failure-rate interval [$COARSE_MIN_FR, $COARSE_MAX_FR]" >&2
  exit 2
fi

{
  echo "coarse_result=found"
  echo "coarse_start_width=$FOUND_W"
  echo "coarse_start_failure_rate=$FOUND_FR"
} >> "$AUTO_SUMMARY"

STOP_W=""
STOP_FR=""
for (( w=FOUND_W; w<=MAX_W; ++w )); do
  cmd=(
    "$RUN_BIN"
    -u 0
    -nn "$NN"
    -ratio "$RATIO"
    -logp "$LOGP"
    -w0 "$w"
    -w1 "$((w+1))"
    -th_cnt "$TH_CNT"
    -trials "$FULL_TRIALS"
  )
  if ((${#EXTRA_ARGS[@]})); then
    cmd+=("${EXTRA_ARGS[@]}")
  fi

  result="$(run_one full "$w" "${cmd[@]}")"
  IFS=$'\t' read -r fr log2_fr failures trials avg_ms <<< "$result"

  status="continue"
  if float_lt "$fr" "$TARGET_FR"; then
    status="stop"
    STOP_W="$w"
    STOP_FR="$fr"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$w" "$fr" "$log2_fr" "$failures" "$trials" "$avg_ms" "$status" \
    >> "$AUTO_FULL_TSV"

  echo "[full] w=$w, failure_rate=$(printf '%.6g' "$fr"), ~$(format_rate_pow2 "$fr")"

  if [[ "$status" == "stop" ]]; then
    echo "[full] stop at w=$w, threshold=$(format_rate_pow2 "$TARGET_FR"), observed=$(format_rate_pow2 "$fr")"
    break
  fi
done

WIDTH_TEST_TSV="tuning/okvs_width/${NN}_${ratio_fixed}_${LOGP}_${FULL_TRIALS}.tsv"

if [[ -n "$STOP_W" ]]; then
  {
    echo "full_result=threshold_reached"
    echo "full_stop_width=$STOP_W"
    echo "full_stop_failure_rate=$STOP_FR"
  } >> "$AUTO_SUMMARY"
else
  {
    echo "full_result=max_w_reached"
    echo "reason=target threshold not reached before max_w"
  } >> "$AUTO_SUMMARY"
fi

echo "width_test_data_file=$WIDTH_TEST_TSV" >> "$AUTO_SUMMARY"

awk -F '\t' '
  BEGIN {
    n = 0
    sx = sy = sxx = sxy = 0
  }
  NR > 1 && $2 != "NA" && $2+0 > 0 {
    x = $1 + 0
    y = log($2 + 0) / log(2)
    n++
    sx += x
    sy += y
    sxx += x * x
    sxy += x * y
  }
  END {
    den = n * sxx - sx * sx
    if (n >= 2 && den != 0) {
      slope = (n * sxy - sx * sy) / den
      intercept = (sy - slope * sx) / n
      printf "fit_points=%d\nfit_slope=%.17g\nfit_intercept=%.17g\nfit_model=log2_failure_rate ~= %.17g * w + %.17g\n", n, slope, intercept, slope, intercept
    } else {
      printf "fit_points=%d\nfit_slope=NA\nfit_intercept=NA\nfit_model=NA\n", n
    }
  }
' "$AUTO_FULL_TSV" >> "$AUTO_SUMMARY"

FIT_SLOPE="$(awk -F= '/^fit_slope=/{print $2}' "$AUTO_SUMMARY" | tail -n 1)"
FIT_INTERCEPT="$(awk -F= '/^fit_intercept=/{print $2}' "$AUTO_SUMMARY" | tail -n 1)"

echo "[done] autosweep summary: $AUTO_SUMMARY"
echo "[done] full table: $AUTO_FULL_TSV"
echo "[done] width_test TSV: $WIDTH_TEST_TSV"
if [[ -n "$FIT_SLOPE" && "$FIT_SLOPE" != "NA" && -n "$FIT_INTERCEPT" && "$FIT_INTERCEPT" != "NA" ]]; then
  echo "[done] fit: log2_failure_rate ~= ${FIT_SLOPE} * w + ${FIT_INTERCEPT}"
fi
