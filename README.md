### Build and Run

The library can be cloned and built with networking support as
```
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
cd ../..
cmake -S . -B build
cmake --build build
```

This will generate `run` executable in `build` directory. 

To run (DEL-free) PSU with set size n = 2^20, 
```
./run -u 7 -nn 20
```
To see other implementations,
```
./run -list
```
and refer cpp files in the `test` directory.
