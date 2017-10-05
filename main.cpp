#define STORE_TRACKS 
//#define CALCULATE_ENERGY
#define CALCULATE_MOM

#include <iostream>
#include <array>
#include <random>
#include <fstream>
#include <ctime>
#include <math.h> 	// For fabs()
#include <sstream>	// for std::stringstream

#include "Functions.h"	
#include "Constants.h"
#include "threevector.h"

bool DontKeepTrack(threevector Pos, threevector Vel, double MASS, double BMag){
	threevector Bhat(0.0,0.0,1.0);
	threevector AccelDir = (Vel.getunit()^Bhat);
	double AngleWithOrigin = asin(Pos.getunit()*AccelDir);
	double VPerp = sqrt(pow(Vel.getx(),2)+pow(Vel.gety(),2));
	double RhoPerp = (MASS*VPerp)/(echarge*BMag);

	if( sqrt(pow(Pos.getx(),2) + pow(Pos.gety(),2)) > 2.0*RhoPerp + 10.0 ){
		return true;
	}
	return false;
}

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\tShow this help message\n\n"
	<< "\t-r,--radius RADIUS\t\t(m), Specify radius of Dust grain\n\n"
	<< "\t-d,--density DENSITY\t\t(m^-^3), Specify density of Dust grain\n\n"
	<< "\t-p,--potential POTENTIAL\t(arb), Specify the potential of Dust grain normalised to electron temperature\n\n"
	<< "\t-b,--bfield BFIELD\t\t(T), Specify the magnetic field (z direction)\n\n"
	<< "\t-e,--etemp ETEMP\t\t(eV), Specify the temperature of plasma electrons\n\n"
	<< "\t-n,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron density\n\n"
	<< "\t-t,--temp TEMP\t\t\t(eV), Specify the temperature of species of interest\n\n"
	<< "\t-s,--spec_charge SPEC_CHARGE\t(arb), Specify the charge of the species of interest\n\n"
	<< "\t-u,--zmaxdebye ZMAXDEBYE\t(arb), Specify the upper limit of simulation domain as number of distances\n\n"
	<< "\t-l,--zmindebye ZMINDEBYE\t(arb), Specify the lower limit of simulation domain as number of distances\n\n"
	<< std::endl;
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

threevector DebyeHuckelField(threevector Position, double Charge, double Radius, double ElectronDensity, double ElectronTemp, double DebyeLength){
	if(Charge==0.0) return threevector(0.0,0.0,0.0);
	threevector Efield = Charge*(1/(4*PI*epsilon0*Position.mag3()))*exp(-(Position.mag3()*Radius)/DebyeLength)
				*(1/Position.mag3()+Radius/DebyeLength)*Position.getunit();

	return Efield;
}

int main(int argc, char* argv[]){
	
	// ***** TIMER AND FILE DECLERATIONS ***** //
	clock_t begin = clock();
	std::string filename = "Data/MagPlasData";
	DECLARE_TRACK();
	std::ofstream AngularMomentumDataFile;	// Data file for angular momentum information


	// ***** DEFINE DUST PARAMETERS ***** //
	double Radius 		= 1e-8;//20e-6;		// m, Radius of dust
	double Density 		= 1850;//19600;//200;			// Be, W, Glass, kg m^(-3), for Tungsten
	double Potential	= -2.5;			// Coulombs, Charge in 
	double BMag 		= 1e-1;//1e-2;			// Tesla, Magnitude of magnetic field

	// ***** DEFINE PLASMA PARAMETERS ***** //
	double eTemp 	= 1.0;	// Electron Temperature, eV
	double eDensity	= 1e18;//1e14;	// m^(-3), Electron density
	double TEMP 	= 1.0;	// eV, This is the Temperature of the species being considered
	int SPEC_CHARGE	= 1.0;	// arb, This is the charge of the species, should be +1.0 or -1.0 normally 
	double MASS;		// kg, This is the mass of the Species being considered
	double zMaxDebye	= 6.0;			// Arb, Number of debye distances max of simulation is
	double zMinDebye	= 6.0;			// Arb, Number of debye distances max of simulation is


	// ***** DETERMINE USER INPUT ***** //
	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 	|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--radius" 	|| arg == "-r" ) 	InputFunction(argc,argv,i,ss0,Radius);
		else if( arg == "--density" 	|| arg == "-d" )	InputFunction(argc,argv,i,ss0,Density);
		else if( arg == "--potential" 	|| arg == "-p" )	InputFunction(argc,argv,i,ss0,Potential);
		else if( arg == "--bfield" 	|| arg == "-b" )	InputFunction(argc,argv,i,ss0,BMag);
		else if( arg == "--etemp" 	|| arg == "-e" )	InputFunction(argc,argv,i,ss0,eTemp);
		else if( arg == "--edensity" 	|| arg == "-n" )	InputFunction(argc,argv,i,ss0,eDensity);
		else if( arg == "--temp" 	|| arg == "-t" )	InputFunction(argc,argv,i,ss0,TEMP);
		else if( arg == "--spec_charge" || arg == "-s" )	InputFunction(argc,argv,i,ss0,SPEC_CHARGE);
		else if( arg == "--zmaxdebye" 	|| arg == "-u" )	InputFunction(argc,argv,i,ss0,zMaxDebye);
		else if( arg == "--zmindebye" 	|| arg == "-l" )	InputFunction(argc,argv,i,ss0,zMinDebye);
                else{
			sources.push_back(argv[i]);
		}
	}

	// If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
	if( SPEC_CHARGE > 0 )	MASS 	= Mp;
	else{
		MASS	= Me;
		eTemp	= TEMP;
	}

	// ***** DEFINE FIELD PARAMETERS ***** //
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.
	threevector BField = BMag*Bhat;
	double DebyeLength 	= sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2)));			// m, Debye Length
//	double Charge 		= SPEC_CHARGE*eTemp*Potential*(4*PI*epsilon0*Radius)*exp(-Radius/DebyeLength); 	// Coulombs, Charge
	double Charge 		= SPEC_CHARGE*eTemp*Potential*(4*PI*epsilon0*Radius); 	// Coulombs, Charge
	double StandardDev=(TEMP*echarge)/(MASS*Radius);	// Thermal Velocity normalised to Dust Grain radii
	double RhoTherm = (MASS*sqrt(StandardDev))/(echarge*BMag);	// Thermal GyroRadius normalised to dust grain radii


	// ***** DEFINE SIMULATION SPACE ***** //
//	double zmax 	= 1.0+zMaxDebye*DebyeLength/Radius;	// Top of Simulation Domain, in Dust Radii
//	double zmin 	= -1.0-zMinDebye*DebyeLength/Radius;	// Bottom of Simulation Domain, in Dust Radii
	double DebyeHuckelImpactParameter = (DebyeLength/Radius)*LambertW((Radius*pow(naturale,2))/DebyeLength);
//	double CoulombImpactParameter	= sqrt(fabs(Charge/(4*PI*epsilon0*BMag)))/Radius;
//	double CoulombImpactParameter	= sqrt(fabs(Charge/(4*PI*epsilon0*StandardDev*Radius)))*pow(epsilon0/MASS,0.25)/Radius;
	double CoulombImpactParameter	= 24.0;
	double zmax 	= 1.0+CoulombImpactParameter;	// Top of Simulation Domain, in Dust Radii
	double zmin 	= -1.0-CoulombImpactParameter;	// Bottom of Simulation Domain, in Dust Radii
	std::cout << "\nRhoTherm = " << RhoTherm;
	std::cout << "\nDHIP = " << DebyeHuckelImpactParameter;
	std::cout << "\nCIP = " << CoulombImpactParameter; std::cin.get();
	double TimeStep = 1e-12;
	double TimeStepTemp = 2*PI*MASS/(echarge*BMag*1000);
	TimeStep = TimeStepTemp;


	// ***** DEFINE RANDOM NUMBER GENERATOR ***** //
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
					<<Radius<<"\t"<<Density<<"\t"<<Charge<<"\t"<<BMag << "\t" << DebyeLength/Radius << "\n";
	AngularMomentumDataFile << "#Output:\tIonNumber\tCollectedIons\tAngularVelocity\tAngularMomentum\n\n";


	// ***** BEGIN LOOP OVER PARTICLE ORBITS ***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	threevector TotalAngularMom(0.0,0.0,0.0);
	DECLARE_LMSUM();
	DECLARE_AMSUM();
	unsigned int j(0), i(0);
	for( i; i < 1000; i ++){

		// ***** RANDOMISE VELOCITY ***** //
		// If the parallel velocity is < vmin, we lose energy which leads to orbits deviating. So we don't include these
		std::normal_distribution<double> Gaussdist(0,sqrt(StandardDev));
		threevector Velocity(Gaussdist(mt),Gaussdist(mt),-abs(Gaussdist(mt)));	// Start with negative z-velocity
		
		INITIAL_VEL();		// For energy calculations

		threevector MeanPosition(0.0,0.0,0.0);
		threevector MeanVelocity(0.0,0.0,0.0);
		threevector OldVelocity(0.0,0.0,0.0);
		threevector OldPosition(0.0,0.0,0.0);

		// ***** RANDOMISE POSITION ***** //
		double ImpactParameter = (1.0+2.0*RhoTherm+CoulombImpactParameter);
		std::uniform_real_distribution<double> dist(-ImpactParameter, ImpactParameter); // IONS
		threevector Position(dist(mt),dist(mt),zmax);
		if( DontKeepTrack(Position,Velocity,MASS,BMag) )
			continue;

		OPEN_TRACK(filename + std::to_string(i) + ".txt");

		// ***** DO PARTICLE PATH INTEGRATION ***** //

		// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
//		threevector EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength)*(1/pow(Radius,2));
		threevector EField = CoulombField(Position,Charge)*(1/pow(Radius,2));
		RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);

//		std::cout << "\n" << Position << "\t" << Velocity;
		UpdateVelocityBoris(MASS,EField,BField,-0.5*TimeStep,Velocity);

		// While the particle is not inside the sphere, not outside the simulation domain
		// And not over a specified number of iterations to catch trapped orbits
		int p(0);
		while( Position.mag3() > 1.0  && Position.getz() > zmin && Position.getz() <= zmax ){
//			EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength)*(1/pow(Radius,2));
			EField = CoulombField(Position,Charge)*(1/pow(Radius,2));
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
		}else{
			// Calculate momentum changes to determine momentum given to dust grain
			FINAL_VEL(); 
			OUTPUT_VEL(i); OUTPUT_VEL("\t"); 
			OUTPUT_VEL(100*(FinalVelMag.mag3()/InitialVelMag.mag3()-1));  OUTPUT_VEL("\n");
			
			threevector CylindricalRadius(Position.getx(),Position.gety(),0.0);
			LinearMomentumSum += MASS*(FinalVelMag-InitialVelMag);
			AngularMomentumSum += MASS*(CylindricalRadius^(FinalVelMag-InitialVelMag));
//			OUTPUT_MOM(i); OUTPUT_MOM("\t");
//			OUTPUT_MOM( MASS*(FinalVelMag-InitialVelMag) ); OUTPUT_MOM("\t");
//			OUTPUT_MOM( MASS*(CylindricalRadius^(FinalVelMag-InitialVelMag)) ); OUTPUT_MOM("\t");
//			OUTPUT_MOM( LinearMomentumSum*(1.0/i) ); OUTPUT_MOM("\t");
//			OUTPUT_MOM( AngularMomentumSum*(1.0/i) ); OUTPUT_MOM("\n");
		}

		CLOSE_TRACK();
//		std::cout << "\ni :\t" << i;
	}
	AngularMomentumDataFile << "\n\n" << i << "\t" << LinearMomentumSum << "\t" << AngularMomentumSum;
//	OUTPUT_MOM(i); OUTPUT_MOM("\t");
//	OUTPUT_MOM( LinearMomentumSum ); OUTPUT_MOM("\t");
//	OUTPUT_MOM( LinearMomentumSum*(1.0/i) ); OUTPUT_MOM("\t");
//	OUTPUT_MOM( AngularMomentumSum ); OUTPUT_MOM("\t");
//	OUTPUT_MOM( AngularMomentumSum*(1.0/i) ); OUTPUT_MOM("\n");
//	std::cout << "\n\nTotalAngularVel = " << TotalAngularVel;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

	AngularMomentumDataFile.close();
//	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	return 0;
}
