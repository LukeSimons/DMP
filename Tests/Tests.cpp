#include <ctime>			// Time program
#include <sstream>			// std::stringstream
#include <string>			// std::string
#include <vector>			// std::vector
#include <random>    			//!< for std::normal_distribution<> etc.


#include "rand_mwts.h"
#include "MotionTest.h" // CHARGING TESTS

using namespace dimplconsts;

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\tShow this help message\n\n"
	<< "\t-m,--mode MODE\tstring variable defining the Test to be run. Available Options:\n"
	<< "\t\tEField			   : Test motion of charge in central electric field, comparing to Euler method\n"
	<< "\t\tVelocityDistribution	   : Test if velocity distribution is maxwellian at the centre of the simulation\n"
	<< "\t-c,--charge CHARGE\tdouble variable defining the number of charges on the dust grain\n"
	<< "\t-t,--timestep TIMESTEP\tdouble variable defining the normalised time step of the simulation\n"
	<< "\t-f,--filename FILENAME\tstring variable defining the filename to be compared to in the test:\n";
}

template<typename T> int InputFunction(int &argc, char* argv[], int &i, std::stringstream &ss0, T &Temp){
	if (i + 1 < argc) { // Make sure we aren't at the end of argv!
		i+=1;
		ss0 << argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
		ss0 >> Temp;
		ss0.clear(); ss0.str("");
		return 0;
	}else{ // Uh-oh, there was no argument to the destination option.
		std::cerr << "\noption requires argument." << std::endl;
		return 1;
	}
}

int main(int argc, char* argv[]){

	// Determine user input for testing mode
	std::string Test_Mode("");
	std::string File_Name("");
	std::vector <std::string> sources;
	std::stringstream ss0;
	double Charge (0.0);
	double Time_Step(0.01);
	int Seed(1);

	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 	|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--mode"	|| arg == "-m" )	InputFunction(argc,argv,i,ss0,Test_Mode);
		else if( arg == "--charge"	|| arg == "-c" )	InputFunction(argc,argv,i,ss0,Charge);
		else if( arg == "--timestep"	|| arg == "-t" )	InputFunction(argc,argv,i,ss0,Time_Step);
		else if( arg == "--filename"	|| arg == "-f" )	InputFunction(argc,argv,i,ss0,File_Name);
		else if( arg == "--seed"	|| arg == "-se" )	InputFunction(argc,argv,i,ss0,Seed);
		else{
			sources.push_back(argv[i]);
		}
	}

	std::cout.precision(10);
	int out(0);
	if( Test_Mode == "EField" ){
		// Test bout 1,
		// Test if motion of particle in constant central electric field in DiMPle
		// is the same as a numerical solution for charged particle motion in central electric field with Euler method.
		if( Charge == 0.0 || Time_Step == 0.0 ){
			std::cout << "Warning! Charge = 0 or Time_Step =0, -c and -t are required inputs for this test!\n";
		}
		out = MotionTest(Test_Mode,File_Name,Time_Step,Charge);
	}else if( Test_Mode == "VelocityDistribution" ){
		// Test bout 2,
		// Test if velocity distribution is maxwellian at the centre of the simulation

		if( Time_Step == 0.0 ){
			std::cout << "Warning! Charge = 0 or Time_Step =0, -c and -t are required inputs for this test!\n";
		}
		
		std::vector<double> zpos, vz;
		double zmin=-10;
		double zmax=10;
		double DriftNorm=0;
		double ThermalVel=1;
		std::mt19937 mt(Seed);
		double ProbUpper = 0.5;
		double TotalTime = 0.0;
		double MaxTime = 200;
    		std::uniform_real_distribution<double> rad(0.0, 1.0);
        	std::normal_distribution<double> Gaussdist(0.0,ThermalVel);
		while( TotalTime < MaxTime ){
			for( int i(0); i < 10; i ++ ){

				// Generate z-velocity from Maxwellian Distribution
				vz.push_back(Gaussdist(mt)); 		
		        	if( rad(mt) > ProbUpper ){
					zpos.push_back(zmax);

					// Generate z-velocity from Makkonen one-way flux
					// vz.push_back(-rand_mwts(-DriftNorm,ThermalVel,mt));
					if( vz[i] > 0 )
						vz[i] = -vz[i];


			        }else{
					zpos.push_back(zmin);
					// Generate z-velocity from Makkonen one-way flux
					// vz.push_back(rand_mwts(DriftNorm,ThermalVel,mt)); 

					// Generate z-velocity from Maxwellian Distribution
					if( vz[i] < 0 )
						vz[i] = -vz[i];
			        }
				
			}
			std::cout << "\nTime = " << TotalTime << "\n";
			TotalTime += Time_Step;
			for( unsigned int j(0); j < zpos.size(); j ++){
				zpos[j] += vz[j]*Time_Step;
				if( zpos[j] > zmax || zpos[j] < zmin ){
					zpos[j] = 0;
					vz[j] = 0;
				}else if( zpos[j] < 0.5 && zpos[j] > -0.5 ){
					if( zpos[j] != 0 && TotalTime > 199)
						std::cout << vz[j] << ", " << zpos[j] << "\n";
				}

			}
			std::cout << "\n";
		}


	}else{
		std::cout << "\n\nInput not recognised! Exiting program\n.";
	}

	if( out == 1 ) std::cout << "\n# PASSED!";
	if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
	if( out == -1 ) std::cout << "\n# FAILED!";

	std::cout << "\n\n*****\n\n\n";	
	return 0;
}
