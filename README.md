# RLWE-OKVS

An implementation of the ssPMT protocol from RLWE-batched homomorphic OKVS
decoding, together with the private set operation (PSO) protocols derived
from it: PSI-Cardinality, PSI-Threshold, PSI-Sum, and PSU.

### Getting the source

This anonymized repository does not support `git clone`. Instead, download
the repository archive from the landing page and extract it. The
`thirdparty/` dependencies are git submodules and are not included in the
archive; fetch them at the pinned versions:

```bash
cd thirdparty
git clone https://github.com/microsoft/SEAL.git SEAL
git -C SEAL checkout v4.4.3
git clone https://github.com/osu-crypto/libOTe.git libOTe
git -C libOTe checkout d5c01fcb6afdec6a01f7e7a1433ed8d86868c85a
git -C libOTe submodule update --init --recursive
git clone https://github.com/Cyan4973/xxHash.git xxHash
git -C xxHash checkout 82cead715cbfddd9e6093db8df95155077ce6e64
cd ..
```

### Build

```bash
cd thirdparty/SEAL
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
