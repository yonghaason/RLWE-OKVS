#!/bin/bash
# One-shot WAN benchmark. Creates an isolated pair of network namespaces
# joined by a veth pair (per-direction htb rate + netem delay), runs the
# chosen protocol with one party per namespace, prints both reports, and
# removes the namespaces again. Nothing outside the namespaces is touched.
#
#   sudo ./benchmark.sh [experiment idx] -nn [logn] -bw [bandwidth] [-v] [extra args]
#
#   idx : 1 SSPMT   2 PSU   3 PSI-Card   4 PSI-Sum   5 PSI-Threshold
#   -nn : log2 of the set size (passed through to ./run)
#   -bw : per-direction bandwidth, e.g. 100mbit / 1gbit / 10gbit (default 10gbit)
#           > 1gbit  ->  latency 0.1 ms (netem delay 0.05 ms per direction)
#          <= 1gbit  ->  latency  80 ms (netem delay   40 ms per direction)
#   -v  : verbose -- print parameters, timings, and communication
#
# Example:  sudo ./benchmark.sh 2 -nn 20 -bw 100mbit -v

set -e

NS0=kkls0; NS1=kkls1
V0=kklsv0; V1=kklsv1
IP0=10.99.0.1; IP1=10.99.0.2
PORT=1212
LIMIT=5000

SCRIPT_DIR=$(dirname "$(realpath "$0")")
BIN="$SCRIPT_DIR/build/run"
LOGDIR="$SCRIPT_DIR/log/benchmark"
RUNUSER=${SUDO_USER:-root}

usage() { grep "^#" "$0" | head -17 | sed 's/^# \?//'; exit 1; }

if [ "$(id -u)" -ne 0 ]; then echo "please run with sudo"; usage; fi
[ -x "$BIN" ] || { echo "binary not found: $BIN (build first)"; exit 1; }

IDX="${1:-}"
[[ "$IDX" =~ ^[1-5]$ ]] || usage
shift

# Pull -bw out of the arguments; everything else goes to ./run verbatim.
BW=10gbit
ARGS=()
while [ $# -gt 0 ]; do
    if [ "$1" = "-bw" ]; then
        BW="$2"; shift 2
    else
        ARGS+=("$1"); shift
    fi
done

# Normalize the bandwidth to a tc rate string and to bits/s.
bwlc=$(echo "$BW" | tr "A-Z" "a-z")
num=$(echo "$bwlc" | grep -oE "^[0-9.]+")
[ -n "$num" ] || { echo "cannot parse bandwidth: $BW"; exit 1; }
case "$bwlc" in
    *gbit|*g) MULT=1000000000 ;;
    *mbit|*m) MULT=1000000 ;;
    *kbit|*k) MULT=1000 ;;
    *bit|*[0-9]) MULT=1 ;;
    *) echo "cannot parse bandwidth: $BW"; exit 1 ;;
esac
case "$bwlc" in
    *bit) RATE=$bwlc ;;
    *)    RATE="${bwlc}bit" ;;
esac
if [ "$(awk -v n="$num" -v m="$MULT" 'BEGIN{print (n*m > 1e9) ? 1 : 0}')" = 1 ]; then
    DELAY=0.05ms; LAT=0.1ms
else
    DELAY=40ms; LAT=80ms
fi

cleanup() {
    ip netns pids $NS0 2>/dev/null | xargs -r kill 2>/dev/null || true
    ip netns pids $NS1 2>/dev/null | xargs -r kill 2>/dev/null || true
    ip netns del $NS0 2>/dev/null || true
    ip netns del $NS1 2>/dev/null || true
}
trap cleanup EXIT
cleanup   # clear leftovers from a previous crashed invocation

# --- setup ---
ip netns add $NS0
ip netns add $NS1
ip link add $V0 type veth peer name $V1
ip link set $V0 netns $NS0
ip link set $V1 netns $NS1
ip -n $NS0 addr add $IP0/24 dev $V0
ip -n $NS1 addr add $IP1/24 dev $V1
ip -n $NS0 link set lo up
ip -n $NS1 link set lo up
ip -n $NS0 link set $V0 up
ip -n $NS1 link set $V1 up
for pair in "$NS0 $V0" "$NS1 $V1"; do
    set -- $pair
    ip netns exec $1 tc qdisc add dev $2 root handle 1: htb default 10
    ip netns exec $1 tc class add dev $2 parent 1: classid 1:10 htb rate $RATE ceil $RATE
    ip netns exec $1 tc qdisc add dev $2 parent 1:10 handle 10: netem delay $DELAY limit $LIMIT
done
echo "link: $RATE per direction, latency $LAT (RTT $(awk -v d="${DELAY%ms}" 'BEGIN{print 2*d}') ms)"

# --- run ---
mkdir -p "$LOGDIR"; chown $RUNUSER "$LOGDIR" 2>/dev/null || true
SLOG="$LOGDIR/last_sender.log"; RLOG="$LOGDIR/last_recver.log"
set +e
ip netns exec $NS0 sudo -u $RUNUSER "$BIN" -u $IDX \
    -role sender -ip $IP0:$PORT "${ARGS[@]}" > "$SLOG" 2>&1 &
SPID=$!
sleep 0.5
ip netns exec $NS1 sudo -u $RUNUSER "$BIN" -u $IDX \
    -role recver -ip $IP0:$PORT "${ARGS[@]}" > "$RLOG" 2>&1
RRC=$?
wait $SPID
SRC=$?
set -e
chown $RUNUSER "$SLOG" "$RLOG" 2>/dev/null || true
echo "==================== sender ===================="
cat "$SLOG"
echo "==================== recver ===================="
cat "$RLOG"
if [ $SRC -ne 0 ]; then echo "!! sender exited with code $SRC"; fi
if [ $RRC -ne 0 ]; then echo "!! recver exited with code $RRC"; fi

# namespaces are removed by the EXIT trap
