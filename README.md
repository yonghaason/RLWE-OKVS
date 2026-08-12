## This branch: homomorphic-rotation decode (design)

This branch holds the design work for replacing the pre-rotated
(sequenced) BGV decode with receiver-side homomorphic rotations: the
receiver ships two Galois keys (row `+1` and column swap) and the sender
rotates during the encrypted decode, so the shifted copies never travel.
[docs/he-rotation-note.tex](docs/he-rotation-note.tex) Section 3 gives the
dense two-row BGV construction and the rotation-schedule analysis —
input-swap (`m/N + (w-1)` rotations) vs output-swap (`L/2 + (w-1)/2`,
the default), cutting communication by ~38% against the pre-rotated
layout's optimum at `eps = 0.2`. See [TODO.md](TODO.md) for the port plan
and why it is parked (noise budget + parameter re-sweep).

---

### Build

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

### Running Benchmark

| index | protocol |
|---|---|
| 1 | ssPMT |
| 2 | PSU |
| 3 | PSI-Card |
| 4 | PSI-Sum |
| 5 | PSI-Threshold |

```bash
sudo ./benchmark.sh [experiment idx] -nn [logn] -bw [bandwidth] -v
```

This runs the chosen protocol with set size `n = 2^logn`, one party per
process, connected by an emulated WAN link (created automatically and
removed when the run finishes).

- `-bw` sets the per-direction bandwidth, e.g. `100mbit`, `1gbit`, `10gbit`
  (default `10gbit`). Above 1 Gbit the link latency is 0.1 ms; at or below
  1 Gbit it is 80 ms.
- `-v` is verbose: prints the underlying parameters, timings, and
  communication.

For example, PSU with `n = 2^20` over a 100 Mbit WAN:

```bash
sudo ./benchmark.sh 2 -nn 20 -bw 100mbit -v
```
