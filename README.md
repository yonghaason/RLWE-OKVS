### Build and Run

The library can be cloned and built with networking support as follows.

```bash
git clone https://github.com/yonghaason/RLWE-OKVS.git --recursive
cd RLWE-OKVS/thirdparty/SEAL
cmake -S . -B build \
 -DSEAL_USE_INTEL_HEXL=ON \
 -DCMAKE_INSTALL_PREFIX=../install/SEAL
cmake --build build
cmake --install build
cd ../libOTe
python3 build.py --all --boost --sodium -DENABLE_CIRCUITS=ON -DENABLE_LOGVOLE=OFF --install=../install/libOTe
cd ../xxHash
make -j
cd ../..
cmake -S . -B build
cmake --build build
```

This generates the `run` executable in the `build` directory.

### Running Tests (single process)

To see the currently registered implementations, run:

```bash
./run -list
```

You can run each test with `run -u <index>`. Indices 0–6 run both parties
in one process over localhost; for example, PSU with set size `n = 2^20`:

```bash
./run -u 6 -nn 20 -v
```

(`-v` prints the underlying parameters, timings, and communication.)

### WAN Experiments (two processes over an emulated network)

Each protocol also has a `*_net_test` variant that runs ONE party per
process over a real TCP connection, so the link between them can be shaped:

| index | test |
|---|---|
| 7 | SSPMT |
| 8 | PSU |
| 9 | PSI-Card |
| 10 | PSI-Sum |
| 11 | PSI-Threshold |

Both processes regenerate the same fixed-seed inputs, and the receiver
verifies the output, so no extra coordination is needed. Manual usage
(two terminals; the sender doubles as the TCP server and listens on `-ip`):

```bash
./run -u 8 -role sender -ip 127.0.0.1:1212 -nn 20 -v
./run -u 8 -role recver -ip 127.0.0.1:1212 -nn 20 -v
```

For WAN measurements use `experiments/wan_netns.sh`, which creates two
network namespaces joined by a veth pair and shapes **each direction
independently** (htb rate + netem delay) — the same structure as a real
full-duplex WAN link. It never touches the host's interfaces, so other
users of the machine are unaffected, and each `run` flushes the
per-namespace TCP metrics so repetitions are independent.

```bash
sudo experiments/wan_netns.sh setup                # 100mbit, 40ms/direction (RTT 80ms)
sudo experiments/wan_netns.sh setup 100mbit 20ms   # RTT 40ms instead
sudo experiments/wan_netns.sh run 8 -nn 20 -v      # PSU under the emulated WAN
sudo experiments/wan_netns.sh status               # queue/drop statistics
sudo experiments/wan_netns.sh teardown             # remove everything
```

`setup` is needed once; the namespaces persist (at zero cost) until
`teardown` or reboot, and `run` can be repeated freely. To change the
rate/delay, `teardown` and `setup` again. Per-run logs are written to
`log/wan_netns/`.

Do NOT emulate a WAN by putting a single `tc netem`/`htb` on `lo`: both
directions (data and ACKs) then share one FIFO and one rate budget, which
inflates the RTT by an order of magnitude under load and distorts every
timing measurement.
