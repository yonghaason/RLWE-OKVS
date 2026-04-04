#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"
RUN_BIN="${RUN_BIN:-$BUILD_DIR/run}"
PROBE_DIR="${PROBE_DIR:-$ROOT_DIR/tuning}"
RESULTS_DIR="${RESULTS_DIR:-$ROOT_DIR/tuning/realresults}"

NN="${NN:-18}"
TEST_INDEX="${TEST_INDEX:-1}"
OUT_TAG="${OUT_TAG:-real}"
PROBE_FILE="${PROBE_FILE:-}"
PROBE_FILE_WAS_SET=0
if [[ -n "$PROBE_FILE" ]]; then
  PROBE_FILE_WAS_SET=1
fi

declare -a CASES=()

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --nn <value>          Override NN (default: ${NN})
  --test-index <value>  Override test index for PSI_card_test (default: ${TEST_INDEX})
  --out-tag <value>     Output suffix tag (default: ${OUT_TAG})
  --probe-file <path>   Override probe-fit TSV source
  --results-dir <path>  Override output directory
  --help                Show this help

Environment overrides:
  NN, TEST_INDEX, OUT_TAG, BUILD_DIR, RUN_BIN, PROBE_DIR, PROBE_FILE, RESULTS_DIR

Examples:
  $(basename "$0")
  NN=22 $(basename "$0")
  $(basename "$0") --nn 22 --out-tag real
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --nn)
      NN="$2"
      shift 2
      ;;
    --test-index)
      TEST_INDEX="$2"
      shift 2
      ;;
    --out-tag)
      OUT_TAG="$2"
      shift 2
      ;;
    --probe-file)
      PROBE_FILE="$2"
      PROBE_FILE_WAS_SET=1
      shift 2
      ;;
    --results-dir)
      RESULTS_DIR="$2"
      shift 2
      ;;
    --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ "$PROBE_FILE_WAS_SET" -eq 0 ]]; then
  PROBE_FILE="$PROBE_DIR/okvs_nn${NN}_probe_fit.tsv"
fi

if [[ ! -x "$RUN_BIN" ]]; then
  echo "run binary not found or not executable: $RUN_BIN" >&2
  exit 1
fi

if [[ ! -f "$PROBE_FILE" ]]; then
  echo "probe fit file not found: $PROBE_FILE" >&2
  exit 1
fi

while IFS=$'\t' read -r ratio _ _ _ _ _ _ extrapolated _; do
  if [[ "$ratio" == "ratio" || -z "$ratio" ]]; then
    continue
  fi
  rounded_w="$(awk -v x="$extrapolated" 'BEGIN { printf "%d", int(x + 0.5) }')"
  CASES+=("${rounded_w} ${ratio}")
done < <(tail -n +2 "$PROBE_FILE" | sort -t $'\t' -k1,1g)

if [[ ${#CASES[@]} -eq 0 ]]; then
  echo "no cases were loaded from: $PROBE_FILE" >&2
  exit 1
fi

mkdir -p "$RESULTS_DIR"

RAW_FILE="$RESULTS_DIR/nn${NN}_cli_all_raw_${OUT_TAG}.txt"
SUMMARY_FILE="$RESULTS_DIR/nn${NN}_cli_all_summary_${OUT_TAG}.txt"

rm -f "$RAW_FILE" "$SUMMARY_FILE"

for spec in "${CASES[@]}"; do
  set -- $spec
  w="$1"
  m_r="$2"

  printf '=== w=%s m_r=%s test_u=%s nn=%s ===\n' \
    "$w" "$m_r" "$TEST_INDEX" "$NN" >> "$RAW_FILE"
  "$RUN_BIN" -u "$TEST_INDEX" -nn "$NN" -v -w "$w" -m_r "$m_r" >> "$RAW_FILE"
  printf '\n' >> "$RAW_FILE"
done

perl -ne '
BEGIN {
  print join("\t", qw(
    bandWidth bandExpansion wrap comm_mb elapsed_ms
    recv_ctxt_ms decode_send_ms okvs_enc_ms enc_send_ms recv_back_ms decrypt_ms
  )), "\n";
}
if (/^=== w=(\S+) m_r=(\S+) /) {
  ($w, $mr) = ($1, $2);
  ($wrap, $comm, $elapsed, $recv_ctxt, $decsend, $okvs, $encsend, $recvback, $decrypt) = (q{}) x 9;
}
elsif (/^wrap:\s+(\S+)/) { $wrap = $1; }
elsif (/^Sender::Recv ctxts & Serialize\s+([0-9]+\.[0-9]+)/) { $recv_ctxt = $1; }
elsif (/^Sender::Encrypted OKVS Decoding & Send Back\s+([0-9]+\.[0-9]+)/) { $decsend = $1; }
elsif (/^Receiver::OKVS Encoding\s+([0-9]+\.[0-9]+)/) { $okvs = $1; }
elsif (/^Receiver::Encryption & Send\s+([0-9]+\.[0-9]+)/) { $encsend = $1; }
elsif (/^Receiver::Recv back and Serialize\s+([0-9]+\.[0-9]+)/) { $recvback = $1; }
elsif (/^Receiver::Decrypt\s+([0-9]+\.[0-9]+)/) { $decrypt = $1; }
elsif (/^comm .* = ([0-9]+\.[0-9]+)MB/) { $comm = $1; }
elsif (/Passed\e\[0m\s+([0-9]+)ms/) {
  $elapsed = $1;
  print join("\t", $w, $mr, $wrap, $comm, $elapsed,
             $recv_ctxt, $decsend, $okvs, $encsend, $recvback, $decrypt), "\n";
}
' "$RAW_FILE" > "$SUMMARY_FILE"

{
  head -n 1 "$SUMMARY_FILE"
  tail -n +2 "$SUMMARY_FILE" \
    | awk -F"\t" '!seen[$1 FS $2]++' \
    | sort -t $'\''\t'\'' -k2,2g -k1,1n
} > "${SUMMARY_FILE}.sorted"
mv "${SUMMARY_FILE}.sorted" "$SUMMARY_FILE"

echo "probe_file: $PROBE_FILE"
echo "raw: $RAW_FILE"
echo "summary: $SUMMARY_FILE"
cat "$SUMMARY_FILE"
