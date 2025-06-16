#!/bin/bash

source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2023-11-23

git clone -b v00-06_fix https://github.com/jeyserma/k4SimDelphes.git
cd k4SimDelphes
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install

cd ../install
export PATH=$(pwd)/bin:${PATH}
export LD_LIBRARY_PATH=$(pwd)/lib64:${LD_LIBRARY_PATH}

