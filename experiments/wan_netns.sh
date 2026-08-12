#!/bin/bash
# Isolated WAN emulation for two-process experiments (see sspmt_net_test).
#
# Builds two network namespaces joined by a veth pair -- effectively two
# virtual machines and a WAN cable -- and shapes EACH DIRECTION independently
# (htb rate + netem delay), which a single qdisc on lo cannot do. Nothing
# outside the namespaces is touched: other users of the host are unaffected,
# and teardown removes every trace. TCP metrics caches are per-namespace, so
# runs here do not inherit TCP state from the host or from each other after
# a re-setup.
#
#   sudo ./wan_netns.sh setup [RATE] [DELAY] [LIMIT]  # default: 100mbit 40ms 5000
#   sudo ./wan_netns.sh run [TEST] [./run args...]    # e.g. run 8 -nn 20 -v
#   sudo ./wan_netns.sh status
#   sudo ./wan_netns.sh teardown
#
# TEST is the *_net_test index in tests/UnitTests.cpp (default 7):
#   7 SSPMT   8 PSU   9 PSI-Card   10 PSI-Sum   11 PSI-Threshold
#
# NOTE on DELAY: it is applied once per direction, so RTT = 2 x DELAY.
# DELAY=40ms gives RTT 80ms (same as `netem delay 40ms` on lo); use 20ms
# for RTT 40ms.

set -e

NS0=kkls0            # sender side, 10.99.0.1
NS1=kkls1            # receiver side, 10.99.0.2
V0=kklsv0
V1=kklsv1
IP0=10.99.0.1
IP1=10.99.0.2
PORT=1212
TEST_IDX=7           # SSPMT_net_test in tests/UnitTests.cpp

SCRIPT_DIR=$(dirname "$(realpath "$0")")
REPO=$(dirname "$SCRIPT_DIR")
BIN="$REPO/build/run"
LOGDIR="$REPO/log/wan_netns"
RUNUSER=${SUDO_USER:-root}

need_root() {
    if [ "$(id -u)" -ne 0 ]; then
        echo "please run with sudo"; exit 1
    fi
}

case "${1:-}" in
setup)
    need_root
    RATE=${2:-100mbit}; DELAY=${3:-40ms}; LIMIT=${4:-5000}
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
    echo "ready: $NS0($IP0) <-[$RATE per direction, RTT = 2 x $DELAY]-> $NS1($IP1)"
    ;;

run)
    need_root
    shift
    if [[ "${1:-}" =~ ^[0-9]+$ ]]; then TEST_IDX=$1; shift; fi
    if ! ip netns list | grep -q "^$NS0"; then
        echo "namespaces not set up; run: sudo $0 setup"; exit 1
    fi
    [ -x "$BIN" ] || { echo "binary not found: $BIN"; exit 1; }
    mkdir -p "$LOGDIR"; chown $RUNUSER "$LOGDIR" 2>/dev/null || true
    # Fresh TCP state every run, so repetitions are independent.
    ip netns exec $NS0 ip tcp_metrics flush || true
    ip netns exec $NS1 ip tcp_metrics flush || true

    SLOG="$LOGDIR/last_sender.log"; RLOG="$LOGDIR/last_recver.log"
    ip netns exec $NS0 sudo -u $RUNUSER "$BIN" -u $TEST_IDX \
        -role sender -ip $IP0:$PORT "$@" > "$SLOG" 2>&1 &
    SPID=$!
    sleep 0.5
    ip netns exec $NS1 sudo -u $RUNUSER "$BIN" -u $TEST_IDX \
        -role recver -ip $IP0:$PORT "$@" > "$RLOG" 2>&1
    wait $SPID
    chown $RUNUSER "$SLOG" "$RLOG" 2>/dev/null || true
    echo "==================== sender ===================="
    cat "$SLOG"
    echo "==================== recver ===================="
    cat "$RLOG"
    ;;

status)
    ip netns list | grep -E "$NS0|$NS1" || { echo "not set up"; exit 0; }
    for pair in "$NS0 $V0" "$NS1 $V1"; do
        set -- $pair
        echo "--- $1 ($2) ---"
        ip netns exec $1 tc -s qdisc show dev $2
    done
    ;;

teardown)
    need_root
    ip netns del $NS0 2>/dev/null || true
    ip netns del $NS1 2>/dev/null || true
    echo "removed"
    ;;

*)
    grep "^#" "$0" | head -27 | sed 's/^# \?//'
    ;;
esac
