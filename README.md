### This branch: sequencing zero-padding theory

This branch tracks `main` and additionally carries `docs/`, the write-ups
behind the sequencing layer padding: the sender pads its homomorphically
decoded layers up to a public budget, so that the transmitted ciphertext
count — otherwise a function of the sender's set — leaks nothing and can be
simulated from public parameters alone.

- `docs/sequencing-note.tex` — self-contained paper section:
  - *Sequencing as covering.* A layer is a window of ζ consecutive blocks
    with one slot per bin; Hall's condition collapses to block intervals, so
    minimizing the layer count L\* is an interval-covering program with a
    min–max theorem (`thm:minmax`).
  - *The first-fit sequencer is exact* (`thm:greedy`): the left-to-right
    smallest-anchor pass attains L\* on every input — each failure yields an
    interval certificate (blocking lemma) and the certificates account for
    every opened layer. Machine-verified across 5·10⁵ instances; the former
    exact-sequencer fallback is dead code by this theorem and was removed.
  - *The public layer budget* (`thm:budget`): with probability 1−2^−λ over
    the hashing, L\* ≤ L̄ = max_g ⌈(D(g)+1)(b+ζ−1)/(g+ζ−1)⌉, where D(g) is an
    exact binomial quantile at level 2^−λ/(Nb²). Both parties derive L̄ from
    public parameters without interaction; the sender always transmits
    exactly L̄ layers and aborts otherwise (probability ≤ 2^−λ).
- `docs/lemma-layer-budget.tex` — the same result packaged for the paper: a
  main-body lemma bounding the sequencer's output directly
  (Pr[L > L′] ≤ 2^−λ) plus the appendix proof, with first-fit optimality as
  a remark.

The implementation lives on `main` (`resolveLayerBudget`,
`resolveSpanBlocks`, `sequenceLayers` in `rlwe-okvs/rlwe-okvs/sspmt.{h,cpp}`).

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
