#!/bin/sh

module load intel-suite

rm DiMPl
icc -O3 -std=c++14 -qopenmp -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -o DiMPl main.cpp threevector.cpp Constants.cpp rand_mwts.cpp Field_Point.cpp Field_Map.cpp -I.
