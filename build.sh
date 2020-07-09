
g++ -std=c++14 -fopenmp -o3 -Wall main.cpp Constants.cpp Functions.cpp threevector.cpp rand_mwts.cpp Field_Point.cpp Field_Map.cpp -I. -o main
#g++ -std=c++14 -o3 -Wall Tests/tmp.cpp threevector.cpp Constants.cpp rand_mwts.cpp Field_Point.cpp Field_Map.cpp -I. -o tmp
g++ -std=c++14 -o3 -Wall Tests/Tests.cpp threevector.cpp Constants.cpp rand_mwts.cpp Field_Point.cpp Field_Map.cpp -I. -I/Tests/ -o test
#g++ -std=c++14 -fopenmp -o3 -Wall main.cpp Constants.cpp Functions.cpp threevector.cpp rand_mwts.cpp Field_Point.cpp Field_Map.cpp -I. -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -lhdf5_hl -lhdf5_cpp -lhdf5_hl_cpp -o main
