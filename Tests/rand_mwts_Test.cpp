#include "Constants.h"	// Define Pre-processor Directives and constants
#include "threevector.h"// for threevector class
#include <math.h> 	// for fabs()
#include <iostream>	// for std::cout

#include <random>	// for std::normal_distribution<> etc.
#include "rand_mwts.h"

#include <math.h>

int main(int argc, char* argv[]){

	double Temperature = 1.0;	// Temperature in eV
	int seed = 1.0;
	double DriftNorm = 0.0;
	std::mt19937 randnumbers(seed);
	double ThermalVel = sqrt(echarge*Temperature/Mp);
	double zvel = rand_mwts(DriftNorm,ThermalVel,randnumbers);       // Generate z-velocity here to determine starting position 
	std::cout << "\n# ThermalVel = " << ThermalVel << "\n";
	for( int i=0; i < 10000; i ++){
		zvel = rand_mwts(DriftNorm,ThermalVel,randnumbers);
		std::cout << "\n" << zvel;
	}
}
