#include <iostream>
#include <fstream>	// for std::fstream
#include <sstream>	// std::stringstream
#include <math.h>	// for cos()
#include "threevector.h" 
#include "Constants.h"

threevector CoulombField(threevector Position, double Charge, double COULOMB_NORM){
	threevector Efield = (COULOMB_NORM*Charge/Position.square())*Position.getunit(); // 0.07, 0.1475

	return Efield;
}

int MotionTest(std::string MotionType, std::string filename, double timestep, double ChargesOnSphere){
	clock_t begin = clock();
	// ********************************************************** //

	std::ofstream TestDataFile;	// Data file for containing the run information
	TestDataFile.open(filename);

	// START NUMERICAL READ IN
	std::fstream file;
	file.open("Data/DiMPl_Track_0.txt");
	if( !file.is_open()){
		std::cout << "Problem opening file.\n";
		return -1;
	}
	std::stringstream ss;
	std::string line;
	int imax = 0;
	std::vector <double> DiMPlPositions, DiMPlVelocities;
	double unused, xpos(0.0), ypos(0.0), zpos(0.0), vx(0.0), vy(0.0), vz(0.0);
	threevector EulerVelocity(0.0,0.0,0.0);
	threevector EulerPosition(0.0,0.0,0.0);
	double MAGNETIC = 100.0;	// Magnetic Field Normalisation

	while( getline(file,line) ){
		imax ++;
		ss << line;
		ss >> xpos >> ypos >> zpos >> vx >> vy >> vz >> unused;
		DiMPlPositions.push_back(zpos);
		DiMPlVelocities.push_back(vz);
		ss.clear();
		ss.str(std::string());
		if( imax == 3 ){
			EulerPosition.setx(xpos*1e-6);
			EulerPosition.sety(ypos*1e-6);
			EulerPosition.setz(zpos*1e-6);
			EulerVelocity.setx(vx/(Mp/(1e-6*echarge*MAGNETIC)));
			EulerVelocity.sety(vy/(Mp/(1e-6*echarge*MAGNETIC)));
			EulerVelocity.setz(vz/(Mp/(1e-6*echarge*MAGNETIC)));
		}
	}

	// END NUMERICAL READ IN

	// START NUMERICAL MODEL COMPARISON
	
	

	double SItimestep = timestep*Mp/(echarge*MAGNETIC);
	double SphereCharge = ChargesOnSphere*echarge;
	
	if( MotionType == "EField" ){

		for( int i = 0; i < imax; i ++){
			threevector Acceleration = (echarge/Mp)*CoulombField(EulerPosition,SphereCharge,1.0/(4.0*PI*epsilon0));

			std::cout << i << "\t" << DiMPlPositions[i] << "\t" << DiMPlVelocities[i]
			<< "\t" << EulerPosition*(1.0/1e-6) << "\t" << EulerVelocity*(Mp/(1e-6*echarge*MAGNETIC)) << "\n";
			EulerPosition+=EulerVelocity*SItimestep;
			EulerVelocity+=Acceleration*SItimestep;

		}
	}
	threevector ModelVelocity(0.0,0.0,DiMPlVelocities[imax-1]);
	EulerVelocity = EulerVelocity*(Mp/(1e-6*echarge*MAGNETIC));
	// END NUMERICAL MODEL COMPARISON
	
	double ReturnVal = 0;

	if( (ModelVelocity - EulerVelocity).mag3() == 0 ) 				ReturnVal = 1;
	else if( (( EulerVelocity.getx() - ModelVelocity.getx() ) / EulerVelocity.getx() ) > 0.05 )	ReturnVal = -1;
	else if( (( EulerVelocity.gety() - ModelVelocity.gety() ) / EulerVelocity.gety() ) > 0.05 )	ReturnVal = -1;
	else if( (( EulerVelocity.getz() - ModelVelocity.getz() ) / EulerVelocity.getz() ) > 0.05 )	ReturnVal = -1;
	else										ReturnVal = 2;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelVel = " << ModelVelocity << "m s^-1 : EulerVel = " << EulerVelocity << "m s^-1";

	if( EulerVelocity.mag3() != 0.0 )
		std::cout << "\nPercentage Deviation: Mag = " << 100*fabs(1.0-ModelVelocity.mag3()/EulerVelocity.mag3()) <<"%\n*****\n\n";
	if( EulerVelocity.getx() != 0.0 )
		std::cout << "\nPercentage Deviation: Xdir = " << 100*fabs(1.0-ModelVelocity.getx()/EulerVelocity.getx()) <<"%\n*****\n\n";
	if( EulerVelocity.gety() != 0.0 )
	std::cout << "\nPercentage Deviation: Ydir = " << 100*fabs(1.0-ModelVelocity.gety()/EulerVelocity.gety()) <<"%\n*****\n\n";
	if( EulerVelocity.getz() != 0.0 )
		std::cout << "\nPercentage Deviation: Zdir = " << 100*fabs(1.0-ModelVelocity.getz()/EulerVelocity.getz()) <<"%\n*****\n\n";

	return ReturnVal;
}

