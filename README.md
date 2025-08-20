### Build

The library can be cloned and built with networking support as
```
git clone https://github.com/yonghaason/RLWE-OKVS.git --recursive
cd rlwe-okvs/thirdparty/SEAL
cmake -S . -B build \
 -DSEAL_USE_INTEL_HEXL=ON \
 -DCMAKE_INSTALL_PREFIX=../install/SEAL \
 -DSEAL_BUILD_EXAMPLES=ON \ 
 -DSEAL_BUILD_BENCH=ON
cmake --build build
cmake --install build
cd ../libOTe
python3 build.py --all --boost --sodium --install=../install/libOTe
cd ../..
cmake -S . -B build
cmake --build build
```