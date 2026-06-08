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

### Running Tests

To see the currently registered implementations, run:

```bash
./run -list
```

You can run each test by running `run -u <index>`. 

For example, to run RPMT-based PSU with set size `n = 2^20`:

```bash
./run -u 2 -nn 20 -v
```

(-v is for printing underlying parameters, timings, and communications.)
