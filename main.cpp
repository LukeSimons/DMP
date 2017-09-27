#include <iostream>
#include <array>
#include <random>
#include <fstream>

#include "Constants.h"
#include "threevector.h"

threevector CoulombField(threevector Position){
	double Charge=2000*echarge;//2000*echarge;
	threevector Efield = (Charge/(4*PI*epsilon0*Position.square()))*Position.getunit();
	return Efield;
}

threevector DebyeHuckelField(threevector Position, double Radius){
	double Charge=2000*echarge;
	if(Charge==0.0) return threevector(0.0,0.0,0.0);
	double ElectronTemp=1;		// eV
	double ElectronDensity=1e18;	// m^(-3)
	double DebyeLength = sqrt((epsilon0*Charge*ElectronTemp)/(ElectronDensity*pow(echarge,2)));
	threevector Efield = Charge*(1/(4*PI*epsilon0*Position.mag3()))*exp(-(Position.mag3()*Radius)/DebyeLength )*(1/Position.mag3()+Radius/DebyeLength)*Position.getunit();

	return Efield;
}

int main(){
	clock_t begin = clock();

	std::ofstream DataFile;		// Output data file
	// ***** DEFINE DUST PARAMETERS ***** //
	// Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
/*
	char EmissivityModel = 'c'; 	// Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
	char ExpansionModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
													// and (z)ero expansion
	char HeatCapacityModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
	char BoilingModel = 'y'; 	// Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
	double DustRadius=1e-6;
	double DustTemperature=300;
	std::array<char,4> ConstModels  = {EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	Matter* Dust = new Tungsten(DustRadius,DustTemperature,ConstModels);
	Dust->update_motion(Position,Velocity,Rotation);
*/

	// ***** DEFINE SIMULATIO SPACE ***** //
	double zmax = 10;
	
	// ***** DEFINE DUST PARAMETERS ***** //
	double Radius = 100e-6;
	double Charge=2000*echarge;
	double ElectronTemp=1;		// eV
	double IonTemp=1;		// eV
	double ElectronDensity=1e18;	// m^(-3)
	double DebyeLength = sqrt((epsilon0*Charge*ElectronTemp)/(ElectronDensity*pow(echarge,2)));

	// ***** DEFINE FIELD PARAMETERS ***** //
	double BMag = 1e-2;
	threevector Bhat(0.0,0.0,1.0);

	// ***** DEFINE RANDOM INITIAL CONDITIONS ***** //
	std::random_device rd;
	std::mt19937 mt(rd());


	std::string filename = "Data/MagPlasData";
	for( unsigned int i(0); i < 100; i ++){

		// VELOCITY
		DataFile.open(filename + std::to_string(i) + ".txt");
		double StandardDeve=(ElectronTemp*echarge)/Me;	// FOR ELECTRONS
		double StandardDevi=(IonTemp*echarge)/Mp;	// FOR IONS
		std::normal_distribution<double> Gaussdist(0.0,sqrt(StandardDeve));
		threevector Velocity(Gaussdist(mt),Gaussdist(mt),-abs(Gaussdist(mt)));	// Start with negative z-velocity
		Velocity = Velocity*(1/Radius);			// Normalise speed to dust grain radius
		threevector OldVelocity(0.0,0.0,0.0);

		// POSITION
		double Vperp = sqrt(pow(Velocity.getx(),2)+pow(Velocity.gety(),2));
		double RhoPerp = Me*Vperp/(echarge*BMag);// FOR ELECTRONS
		std::uniform_real_distribution<double> dist(-(1.0+RhoPerp), 1.0+RhoPerp);
		threevector Position(dist(mt),dist(mt),zmax);	// Distance normalised to dust grain size, Start 10 Radii above dust


		double Rotation = 0;

		double TimeStep = 1e-15;
		double TimeStepTemp = 2*PI*Me/(echarge*BMag*10000);
//		std::cout << "TimeStepTemp = " << TimeStepTemp; std::cin.get();
		TimeStep = TimeStepTemp;

		// ***** DO INTEGRATION ***** //

//		std::cout << "\n" << Position << "\t\t" << Velocity; // << "\t\t" << CoulombField(MeanPosition)*(1/pow(Radius,3)) << "\t\t" << (BMag*MeanVelocity^Bhat) <<  "\n";
		DataFile << "\n" << Position << "\t" << Velocity;

		threevector OldPosition = Position;
		Position+=TimeStep*Velocity;
		threevector MeanPosition = 0.5*(Position+OldPosition);
		threevector MeanVelocity = ((Position-OldPosition)*(1/TimeStep));
		threevector DeltaV=TimeStep*(echarge/Me)*(DebyeHuckelField(MeanPosition,Radius)*(1/pow(Radius,3))+BMag*MeanVelocity^Bhat); // FOR ELECTRONS

		Velocity += DeltaV;

		DataFile << "\n" << Position << "\t" << Velocity;

		while(Position.mag3() > 1.0 && Position.getz() > 0.0 && Position.getz() < zmax ){
			OldPosition = Position;
			Position+=TimeStep*Velocity;
			MeanPosition = 0.5*(Position+OldPosition);
			MeanVelocity = ((Position-OldPosition)*(1/TimeStep));
//			DeltaV=(TimeStep*echarge/Me)*(DebyeHuckelField(MeanPosition,Radius)*(1/pow(Radius,3))+ BMag*Velocity.mag3()*(MeanVelocity.getunit()^Bhat)); // FOR ELECTRONS
			DeltaV=(TimeStep*echarge/Me)*(CoulombField(MeanPosition)*(1/pow(Radius,3))+ BMag*Velocity.mag3()*(MeanVelocity.getunit()^Bhat)); // FOR ELECTRONS
		//	std::cout << "\nCoeff : " << TimeStep*(echarge/Me)*BMag;
		//	std::cout << "\nVxB : " << (MeanVelocity^Bhat);


			OldVelocity = Velocity;
			Velocity += DeltaV;
//			threevector DeltaV=(DebyeHuckelField(MeanPosition,Radius)*(1/pow(Radius,3))+BMag*MeanVelocity^Bhat)*(echarge/Mp); // FOR IONS
//			threevector DeltaV=(CoulombField(MeanPosition)*(1/pow(Radius,3))+BMag*MeanVelocity^Bhat)*(echarge/Me);

			DataFile << "\n" << Position << "\t" << Velocity;// << "\t\t" << (TimeStep*(echarge/Me)*BMag*MeanVelocity^Bhat); // << "\t\t" << CoulombField(MeanPosition)*(1/pow(Radius,3)) << "\t\t" << (BMag*MeanVelocity^Bhat) <<  "\n";
		}


		if(Position.mag3() < 1.0){ // In this case it was captured!
			threevector FinalVelocity = 0.5*(OldVelocity+Velocity);
			threevector FinalPosition = MeanPosition;
			threevector CylindricalRadius(FinalPosition.getx(),FinalPosition.gety(),0.0);
			double DistanceFromAxis = CylindricalRadius.mag3();
			threevector Torque = Me*DistanceFromAxis*pow(Radius,2)*(CylindricalRadius^FinalVelocity); // FOR ELECTRONS
//			threevector Torque = Mp*DistanceFromAxis*pow(Radius,2)*(CylindricalRadius^FinalVelocity); // FOR IONS
			std::cout << "\nTorque = " << Torque;
		}
		DataFile.close();
	}
//	std::cout << "\nTimeStep = " << TimeStep;
//	std::cout << "\nPosition.mag3() = " << Position.mag3();
//	std::cout << "\nPosition.getz() = " << Position.getz();
//	std::cout << "\nDebyeLength/Radius = " << DebyeLength/Radius;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	return 0;
}
