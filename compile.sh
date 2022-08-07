#!/usr/bin/bash
# g++ -std=c++20 src/cMchecker.cpp -o build/cMchecker

# make build directory if it does not exist
if [[ ! -d build ]]
then
  mkdir build
fi

# switch to build directory, cmake and make
cd build
cmake ..
make
