#!bin/bash
g++ -std=c++14 -fopenmp -o3 -Wall main.cpp Constants.cpp Functions.cpp threevector.cpp -I. -o main 
