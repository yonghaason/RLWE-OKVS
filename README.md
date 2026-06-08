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
sudo apt install libtool
python3 build.py --all --boost --sodium -DENABLE_CIRCUITS=ON --install=../install/libOTe
cd ../xxHash
make -j
cd ../..
cmake -S . -B build
cmake --build build
```

This generates the `run` executable in the `build` directory.

### Available Tests

To see the currently registered implementations, run:

```bash
./run -list
```

### Running Tests

The executable uses `-u <index>` to select an implementation. Since the numeric index can change when tests are added or removed, first run `./run -list` and use the index printed next to the desired test.

To run RPMT-based PSU with set size `n = 2^20`:

```bash
./run -u 2 -nn 20
```

To run secret-shared PMT-based PSU with set size `n = 2^20`:

```bash
./run -u 3 -nn 20
```

Useful options:

```bash
-n <set-size>      # use an explicit set size
-nn <log-set-size> # use set size n = 2^nn
-nt <threads>      # number of worker threads
-v                 # print parameters, timing, and communication
```

For the secret-shared PMT-based PSU test, the following parameter overrides are also supported:

```bash
-w <band-width>
-m_r <band-expansion>
-seq_span <span-blocks>
```

Refer to the C++ files in the `test` directory for implementation details.
