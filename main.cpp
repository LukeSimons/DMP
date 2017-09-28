//#define STORE_TRACKS 

#include <iostream>
#include <array>
#include <random>
#include <fstream>
#include <ctime>

#include "Constants.h"
#include "threevector.h"

threevector CoulombField(threevector Position, double Charge){
	threevector Efield = (Charge/(4*PI*epsilon0*Position.square()))*Position.getunit();
	return Efield;
}

threevector DebyeHuckelField(threevector Position, double Charge, double Radius, double ElectronDensity, double ElectronTemp){
	if(Charge==0.0) return threevector(0.0,0.0,0.0);
	double DebyeLength = sqrt((epsilon0*echarge*ElectronTemp)/(ElectronDensity*pow(echarge,2)));
	threevector Efield = Charge*(1/(4*PI*epsilon0*Position.mag3()))*exp(-(Position.mag3()*Radius)/DebyeLength)
				*(1/Position.mag3()+Radius/DebyeLength)*Position.getunit();

	return Efield;
}

int main(){
	clock_t begin = clock();
	std::string filename = "Data/MagPlasData";
	DECLARE_TRACK();
	std::ofstream AngularMomentumDataFile;	// Data file for angular momentum information

	// ***** DEFINE PLASMA PARAMETERS ***** //
	double eTemp 	= 1.0;	// Electron Temperature, eV
	double eDensity	= 1e18;	// m^(-3), Electron density
	double TEMP 	= 1.0;	// eV, This is the Temperature of the species being considered
	int SPEC_CHARGE	= 1.0;	// arb, This is the charge of the species, should be +1.0 or -1.0 normally 
	double MASS;		// kg, This is the mass of the Species being considered

	// If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
	if( SPEC_CHARGE > 0 )	MASS 	= Mp;
	else			MASS	= Me;
	
	// ***** DEFINE DUST PARAMETERS ***** //
	double Radius 		= 100e-6;		// m, Radius of dust
	double Density 		= 19600;		// kg m^(-3), for Tungsten
	double Potential	= -2.6;			// Coulombs, Charge in 
	double DebyeLength 	= sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2)));
	double Charge 		= eTemp*Potential*(4*PI*epsilon0*Radius);//*exp(-Radius/DebyeLength)); // Coulombs, Charge in 
//	std::cout << "\nCharge = " << Charge << "C\nElectrons = " << Charge/echarge;
//	std::cout << "e-\nDebyeLength = " << DebyeLength/Radius << "arb\nRadius = " << Radius;
//	std::cout << "m\nPotential = " << Potential << "arb\neTemp = " << eTemp << "ev\n"; std::cin.get();

	// ***** DEFINE SIMULATION SPACE ***** //
	double zmax 	= 1.0+10.0*DebyeLength/Radius;	// Top of Simulation Domain, in Dust Radii
	double zmin 	= -1.0-5.0*DebyeLength/Radius;	// Bottom of Simulation Domain, in Dust Radii

	// ***** DEFINE FIELD PARAMETERS ***** //
	double BMag = 1e-1;		// Tesla, Magnitude of magnetic field
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.

	// ***** DEFINE RANDOM INITIAL CONDITIONS ***** //
	std::random_device rd;		// Create Random Device
	std::mt19937 mt(rd());		// Get Random Method

	// ***** OPEN DATA FILE WITH HEADER ***** //
	time_t now = time(0);		// Get the time of simulation
	char * dt = ctime(&now);
	AngularMomentumDataFile.open(filename + "_AngMom.txt");	
	AngularMomentumDataFile << "## Angular Momentum Data File ##\n";
	AngularMomentumDataFile << "#Date: " << dt;
	AngularMomentumDataFile << "#Input:\tzmax\tzmin\telec_temp\telec_dens\tion_temp\tRadius\tDensity\tCharge\t\tBMag\tDebye\n";
	AngularMomentumDataFile << "#\t"<<zmax<<"\t"<<zmin<<"\t"<<eTemp<<"\t\t"<<eDensity<<"\t\t"<<TEMP<<"\t\t"
					<<Radius<<"\t"<<Density<<"\t"<<Charge<<"\t"<<BMag << "\t" << DebyeLength/Radius << "\n\n";

	// ***** BEGIN LOOP OVER PARTICLE ORBITS ***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	unsigned int j(0);
	for( unsigned int i(0); i < 100; i ++){

		// VELOCITY
		OPEN_TRACK(filename + std::to_string(i) + ".txt");
		double StandardDev=(TEMP*echarge)/MASS;	// FOR IONS
		std::normal_distribution<double> Gaussdist(0.0,sqrt(StandardDev));
		threevector Velocity(Gaussdist(mt),Gaussdist(mt),-abs(Gaussdist(mt)));	// Start with negative z-velocity
		Velocity = Velocity*(1/Radius);			// Normalise speed to dust grain radius
		threevector MeanPosition(0.0,0.0,0.0);
		threevector MeanVelocity(0.0,0.0,0.0);
		threevector OldVelocity(0.0,0.0,0.0);

		// POSITION
		double Vperp = sqrt(pow(Velocity.getx(),2)+pow(Velocity.gety(),2));
		double RhoPerp = (MASS*Vperp)/(echarge*BMag);
		threevector Position(0.0,0.0,zmax);
		if(SPEC_CHARGE > 0.0){
			std::uniform_real_distribution<double> dist(-(1.0+RhoPerp+4*DebyeLength), 1.0+RhoPerp+4*DebyeLength); // FOR IONS
			Position.sety(dist(mt));  // Distance normalised to dust grain size, Start 10 Radii above dust
		}else{
			std::uniform_real_distribution<double> dist(-(1.0+RhoPerp), 1.0+RhoPerp); // FOR ELECTRONS
			Position.sety(dist(mt));  // Distance normalised to dust grain size, Start 10 Radii above dust
		}

		double TimeStep = 1e-12;
		double TimeStepTemp = 2*PI*MASS/(echarge*BMag*10000);
		TimeStep = TimeStepTemp;

		// ***** DO PARTICLE PATH INTEGRATION ***** //

		RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
		// While we're not inside the sphere and we're not outside the simulation domain
		while(Position.mag3() > 1.0 && Position.getz() > zmin && Position.getz() <= zmax ){
			threevector OldPosition = Position;
			Position+=TimeStep*Velocity;
			MeanPosition = 0.5*(Position+OldPosition);
			MeanVelocity = ((Position-OldPosition)*(1/TimeStep));
			// ONLY ONE OF THE FOLLOWING SHOULD BE UNCOMMENTED
//			DebyeHuckel Potential
			threevector DeltaV=((SPEC_CHARGE*TimeStep*echarge)/MASS)
					*(DebyeHuckelField(MeanPosition,Charge,Radius,eDensity,eTemp)*(1/pow(Radius,3))
					+BMag*Velocity.mag3()*(MeanVelocity.getunit()^Bhat));

// 			Coulomb Potential - MAKE SURE TO CHECK POSITION DISTRIBUTION
//			threevector DeltaV=((SPEC_CHARGE*TimeStep*echarge)/MASS)
//					*(CoulombField(MeanPosition,Charge)*(1/pow(Radius,3))
//					+BMag*Velocity.mag3()*(MeanVelocity.getunit()^Bhat));

			OldVelocity = Velocity;
			Velocity += DeltaV;

			RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
//			std::cout << "\n" << Position << "\t" << Velocity;
		}


		if(Position.mag3() < 1.0){ // In this case it was captured!
			threevector FinalVelocity = 0.5*(OldVelocity+Velocity);
			threevector FinalPosition = MeanPosition;
			threevector CylindricalRadius(FinalPosition.getx(),FinalPosition.gety(),0.0);
			double DistanceFromAxis = CylindricalRadius.mag3();
//			threevector AngularMom = MASS*DistanceFromAxis*pow(Radius,2)*(CylindricalRadius^FinalVelocity); 
			// Moment of Inertia for solid sphere
			threevector AngularVel = (15*MASS)/(8*PI*Density*pow(Radius,3))*	
				(CylindricalRadius^(FinalVelocity-AngularVel.mag3()*DistanceFromAxis*AngularVel.getunit()));

			TotalAngularVel += AngularVel;
			j ++;
			std::cout << "\nAngularVel = " << AngularVel << "\nj :\t" << j << "\tTotalAngularVel = " << TotalAngularVel;
			AngularMomentumDataFile << "\nj :\t" << j << "\tTotalAngularVel = " << TotalAngularVel;

		}
		CLOSE_TRACK();
//		std::cout << "\ni :\t" << i;
	}
	std::cout << "\n\nTotalAngularVel = " << TotalAngularVel;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

	AngularMomentumDataFile.close();
	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	return 0;
}
