#!bin/bash
g++ -std=c++14 -fopenmp -o3 -Wall main.cpp Constants.cpp Functions.cpp threevector.cpp rand_mwts.cpp -I. -o main
#g++ -std=c++14 -fopenmp -o3 -Wall main.cpp Constants.cpp Functions.cpp threevector.cpp rand_mwts.cpp -I. -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lhdf5_hl -lhdf5_cpp -lhdf5_hl_cpp -o main
