//#define STORE_TRACKS 

#include <iostream>
#include <array>
#include <random>
#include <fstream>
#include <ctime>

#include "Constants.h"
#include "threevector.h"

/*updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation, p.62*/
static void UpdateVelocityBoris(double MASS, threevector Efield, threevector BField, double dt, threevector &Velocity){

	/*t vector*/
	threevector t = echarge*BField*0.5*dt*(1/MASS);
	
	/*magnitude of t, squared*/
	double t_mag2 = t.square();
	
	/*s vector*/
	threevector s = 2.0*t*(1/(1+t_mag2));
	
	/*v minus*/
	threevector v_minus = Velocity + echarge*Efield*0.5*dt*(1/MASS); 
	
	/*v prime*/
	threevector v_prime = v_minus + (v_minus^t);
	
	/*v prime*/
	threevector v_plus = v_minus + (v_prime^s);
	
	/*v n+1/2*/
	Velocity = v_plus + echarge*Efield*0.5*dt*(1/MASS);
}

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
	double eDensity	= 1e14;	// m^(-3), Electron density
	double TEMP 	= 1.0;	// eV, This is the Temperature of the species being considered
	int SPEC_CHARGE	= 1.0;	// arb, This is the charge of the species, should be +1.0 or -1.0 normally 
	double MASS;		// kg, This is the mass of the Species being considered
	// If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
	if( SPEC_CHARGE > 0 )	MASS 	= Mp;
	else{
		MASS	= Me;
		eTemp	= TEMP;
	}
	
	// ***** DEFINE DUST PARAMETERS ***** //
	double Radius 		= 20e-6;		// m, Radius of dust
	double Density 		= 200;		// kg m^(-3), for Tungsten
	double Potential	= -2.5;			// Coulombs, Charge in 
	double DebyeLength 	= sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2)));
	double Charge 		= eTemp*Potential*(4*PI*epsilon0*Radius)*exp(-Radius/DebyeLength); // Coulombs, Charge in 

	// ***** DEFINE FIELD PARAMETERS ***** //
	double BMag = 1e-2;		// Tesla, Magnitude of magnetic field
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.
	threevector BField = BMag*Bhat;

	// ***** DEFINE SIMULATION SPACE ***** //
	double zmax 	= 1.0+10.0*DebyeLength/Radius;				// Top of Simulation Domain, in Dust Radii
	double zmin 	= -1.0-5.0*DebyeLength/Radius;				// Bottom of Simulation Domain, in Dust Radii
	double vMin	= (Radius*(zmax-zmin)*echarge*BMag)/(MASS*1000);		// THIS MAY NEED TO BE CHANGED FOR DIFFERENT FIELDS

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
	AngularMomentumDataFile << "#Output:\tIonNumber\tCollectedIons\tAngularVelocity\tAngularMomentum\n\n";

	// ***** BEGIN LOOP OVER PARTICLE ORBITS ***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	threevector TotalAngularMom(0.0,0.0,0.0);
	unsigned int j(0);

	for( unsigned int i(0); i < 10; i ++){

		// VELOCITY
		OPEN_TRACK(filename + std::to_string(i) + ".txt");
		double StandardDev=(TEMP*echarge)/MASS;	// FOR IONS
		// If the parallel velocity is < vmin, we lose energy which leads to orbits deviating. So we don't include these
		std::normal_distribution<double> GaussdistXY(0,sqrt(StandardDev));
//		std::normal_distribution<double> GaussdistZ(vMin,sqrt(StandardDev));
		std::normal_distribution<double> GaussdistZ(0,sqrt(StandardDev));
		threevector Velocity(GaussdistXY(mt),GaussdistXY(mt),-abs(GaussdistZ(mt)));	// Start with negative z-velocity
		Velocity = Velocity*(1/Radius);			// Normalise speed to dust grain radius
		
//		double InitialVelMag = Velocity.square();
		threevector MeanPosition(0.0,0.0,0.0);
		threevector MeanVelocity(0.0,0.0,0.0);
		threevector OldVelocity(0.0,0.0,0.0);
		threevector OldPosition(0.0,0.0,0.0);

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
		double TimeStepTemp = 2*PI*MASS/(echarge*BMag*1000);
		TimeStep = TimeStepTemp;

		// ***** DO PARTICLE PATH INTEGRATION ***** //

		// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
		threevector EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp)*(1/pow(Radius,3));
//		threevector EField = CoulombField(Position,Charge)*(1/pow(Radius,3));
		RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);

//		std::cout << "\n" << Position << "\t" << Velocity;
		UpdateVelocityBoris(MASS,EField,BField,-0.5*TimeStep,Velocity);

		// While the particle is not inside the sphere, not outside the simulation domain
		// And not over a specified number of iterations to catch trapped orbits
		int p(0);
		while( Position.mag3() > 1.0  && Position.getz() > zmin && Position.getz() <= zmax ){
// && p < abs(round(5*(zmax-zmin)/(Velocity.getz()*TimeStep)))

			EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp)*(1/pow(Radius,3));
//			EField = CoulombField(Position,Charge)*(1/pow(Radius,3));
			OldVelocity = Velocity;	// For Angular Momentum Calculations
			OldPosition = Position; // For Angular Momentum Calculations
			UpdateVelocityBoris(MASS,EField,BField,TimeStep,Velocity);
			Position+=TimeStep*Velocity;

			if( p%2==0 ){
				RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
			}
			p ++;
//			std::cout << "\n" << Position << "\t" << Velocity;
		}


		if(Position.mag3() < 1.0){ // In this case it was captured!
			threevector FinalVelocity = 0.5*(OldVelocity+Velocity);
			threevector FinalPosition = 0.5*(OldPosition+Position);
			threevector CylindricalRadius(FinalPosition.getx(),FinalPosition.gety(),0.0);
			double DistanceFromAxis = CylindricalRadius.mag3();
			threevector AngularMom = MASS*pow(Radius,2)*(CylindricalRadius^FinalVelocity); 
			// Moment of Inertia for Hollow Sphere
			threevector AngularVel = (9*MASS)/(8*PI*Density*pow(Radius,3))*	
				(CylindricalRadius^(FinalVelocity-AngularVel.mag3()*DistanceFromAxis*AngularVel.getunit()));
			// Moment of Inertia for solid sphere
//			threevector AngularVel = (15*MASS)/(8*PI*Density*pow(Radius,3))*	
//				(CylindricalRadius^(FinalVelocity-AngularVel.mag3()*DistanceFromAxis*AngularVel.getunit()));

			TotalAngularVel += AngularVel;
			TotalAngularMom += AngularMom;
			j ++;
//			std::cout << "\nAngularVel = " << AngularVel << "\nj :\t" << j << "\tTotalAngularVel = " << TotalAngularVel;
			AngularMomentumDataFile << "\n" << i << "\t" << j << "\t" << TotalAngularVel << "\t" << AngularMom;
		}

		CLOSE_TRACK();
//		std::cout << "\ni :\t" << i;
//		double FinalVelMag = Velocity.square();	std::cout << "\nEnergyLoss Percent = " << (FinalVelMag/InitialVelMag-1);
	}
//	std::cout << "\n\nTotalAngularVel = " << TotalAngularVel;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

	AngularMomentumDataFile.close();
//	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	return 0;
}
