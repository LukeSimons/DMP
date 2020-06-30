#!/bin/sh

module load intel-suite
#rm DiMPl
#icc -O3 -std=c++14 -qopenmp -o DiMPl_Debye main.cpp threevector.cpp Constants.cpp rand_mwts.cpp -I.
#icc -O3 -std=c++14 -qopenmp -o DiMPl_Debye main.cpp threevector.cpp Constants.cpp rand_mwts.cpp -I.
icc -O3 -std=c++14 -qopenmp -o DiMPl_Coulomb main.cpp threevector.cpp Constants.cpp rand_mwts.cpp -I.
#icc -O3 -std=c++14 -qopenmp -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2 -o DiMPl_Dynam_Debye main.cpp threevector.cpp Constants.cpp rand_mwts.cpp -I.

#icc -std=c++14 -xCORE-AVX-I -o DMP main.cpp threevector.cpp Constants.cpp -I.
