#!bin/bash
g++ -std=c++14 -fopenmp main.cpp Constants.cpp Functions.cpp threevector.cpp -I. -o main 
