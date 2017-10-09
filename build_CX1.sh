#!/bin/sh

module load intel-suite
icc -O3 -std=c++14 -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -o DMP main.cpp threevector.cpp Constants.cpp -I.

#icc -std=c++14 -xCORE-AVX-I -o DMP main.cpp threevector.cpp Constants.cpp -I.
