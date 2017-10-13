#define STORE_TRACKS 
//#define CALCULATE_ENERGY
#define CALCULATE_MOM

#include <omp.h>	// For parallelisation

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

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\tShow this help message\n\n"
	<< "\t-r,--radius RADIUS\t\t(m), Specify radius of Dust grain\n\n"
	<< "\t-d,--density DENSITY\t\t(m^-^3), Specify density of Dust grain\n\n"
	<< "\t-p,--potential POTENTIAL\t(arb), Specify the potential of Dust grain normalised to electron temperature\n\n"
	<< "\t-m,--magfield MAGFIELD\t\t(T), Specify the magnetic field (z direction)\n\n"
	<< "\t-e,--etemp ETEMP\t\t(eV), Specify the temperature of plasma electrons\n\n"
	<< "\t-n,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron density\n\n"
	<< "\t-t,--temp TEMP\t\t\t(eV), Specify the temperature of species of interest\n\n"
	<< "\t-s,--spec_charge SPEC_CHARGE\t(arb), Specify the charge of the species of interest\n\n"
	<< "\t-u,--zmaxdebye ZMAXDEBYE\t(arb), Specify the upper limit of simulation domain as number of distances\n\n"
	<< "\t-l,--zmindebye ZMINDEBYE\t(arb), Specify the lower limit of simulation domain as number of distances\n\n"
	<< "\t-b,--impactpar IMPACTPAR\t(arb), Specify the radial limit of simulation domain as number of distances\n\n"
	<< "\t-i,--imax IMAX\t(arb), Specify the number of particles to be launched\n\n"
	<< "\t-v,--driftvel DRIFTVEL\t(m s^-^1), Specify the drift velocity of the plasma\n\n"
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
	threevector t = BField*0.5*dt;
	
	/*magnitude of t, squared*/
	double t_mag2 = t.square();
	
	/*s vector*/
	threevector s = 2.0*t*(1/(1+t_mag2));
	
	/*v minus*/
	threevector v_minus = Velocity + Efield*0.5*dt; 
	
	/*v prime*/
	threevector v_prime = v_minus + (v_minus^t);
	
	/*v prime*/
	threevector v_plus = v_minus + (v_prime^s);
	
	/*v n+1/2*/
	Velocity = v_plus + Efield*0.5*dt;
}

threevector CoulombField(threevector Position, double Charge, double e0norm){
	threevector Efield = (Charge/(4*PI*e0norm*Position.square()))*Position.getunit();
	return Efield;
}

threevector DebyeHuckelField(threevector Position, double Charge, double Radius, double ElectronDensity, double ElectronTemp, double DebyeLength, double e0norm){
	if(Charge==0.0) return threevector(0.0,0.0,0.0);
	threevector Efield = Charge*(1.0/(4*PI*e0norm*Position.mag3()))*exp(-Position.mag3()/DebyeLength)
				*(1.0/Position.mag3()+1.0/DebyeLength)*Position.getunit();

	return Efield;
}

int main(int argc, char* argv[]){
	
	// ***** TIMER AND FILE DECLERATIONS ***** //
	clock_t begin = clock();
	std::string filename = "Data/MagPlasData";
	DECLARE_TRACK();
	std::ofstream AngularMomentumDataFile;	// Data file for angular momentum information


	// ***** DEFINE DUST PARAMETERS ***** //
	double Radius 		= 1e-6;	// m, Radius of dust
	double Density 		= 19600;// kg m^-^3, Tungsten
	double Potential	= -2.5;	// Coulombs, Charge in 
	double BMag 		= 1.0; // Tesla, Magnitude of magnetic field


	// ***** DEFINE PLASMA PARAMETERS ***** //
	double eTemp 		= 1.0;	// Electron Temperature, eV
	double eDensity		= 1e18;	//1e14;	// m^(-3), Electron density
	double TEMP 		= 1.0;	// eV, This is the Temperature of the species being considered
	double DriftVel 	= 0.0;	// m s^-1, This is the Temperature of the species being considered
	int SPEC_CHARGE		= 1.0;	// arb, This is the charge of the species, should be +1.0 or -1.0 normally 
	double zMaxDebye	= 50.0;	// Arb, Number of debye distances max of simulation is
	double zMinDebye	= 50.0;	// Arb, Number of debye distances max of simulation is
	double ImpactPar	= 1.0;	// Arb, Multiplicative factor for the Impact Parameter
	unsigned int imax	= 100;

	// ***** DETERMINE USER INPUT ***** //
	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 	|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--radius" 	|| arg == "-r" ) 	InputFunction(argc,argv,i,ss0,Radius);
		else if( arg == "--density" 	|| arg == "-d" )	InputFunction(argc,argv,i,ss0,Density);
		else if( arg == "--potential" 	|| arg == "-p" )	InputFunction(argc,argv,i,ss0,Potential);
		else if( arg == "--magfield" 	|| arg == "-m" )	InputFunction(argc,argv,i,ss0,BMag);
		else if( arg == "--etemp" 	|| arg == "-e" )	InputFunction(argc,argv,i,ss0,eTemp);
		else if( arg == "--edensity" 	|| arg == "-n" )	InputFunction(argc,argv,i,ss0,eDensity);
		else if( arg == "--temp" 	|| arg == "-t" )	InputFunction(argc,argv,i,ss0,TEMP);
		else if( arg == "--spec_charge" || arg == "-s" )	InputFunction(argc,argv,i,ss0,SPEC_CHARGE);
		else if( arg == "--zmaxdebye" 	|| arg == "-u" )	InputFunction(argc,argv,i,ss0,zMaxDebye);
		else if( arg == "--zmindebye" 	|| arg == "-l" )	InputFunction(argc,argv,i,ss0,zMinDebye);
		else if( arg == "--impactpar"	|| arg == "-b" )	InputFunction(argc,argv,i,ss0,ImpactPar);
		else if( arg == "--imax"	|| arg == "-i" )	InputFunction(argc,argv,i,ss0,imax);
		else if( arg == "--driftvel"	|| arg == "-v" )	InputFunction(argc,argv,i,ss0,DriftVel);
                else{
			sources.push_back(argv[i]);
		}
	}

	// If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
	double MASS;			// kg, This is the mass of the Species being considered
	if( SPEC_CHARGE > 0 )	MASS 	= Mp;
	else{
		MASS	= Me;
		eTemp	= TEMP;
	}
	
	// ***** NORMALISATION ***** //
	// Normalise TIME to Gyro period
	// Normalise MASS to Species Mass
	// Normalise DISTANCE to Dust Radius
	// Normalise CHARGE to fundamental charge
        double Tau;
        if( BMag <= 0.0 ){
                Tau     = Radius*sqrt(MASS/(echarge*TEMP));
        }else{
                Tau     = MASS/(echarge*BMag);
        }
	double e0norm 	= epsilon0*MASS*pow(Radius,3)/(pow(echarge*Tau,2));
	double u0norm 	= (4.0*PI*10e-7)/(pow(echarge,2)*MASS*Radius);
	double NormDens = Density*pow(Radius,3)/MASS;
 	double TimeStep	= 1.0/100000;
	double PotNorm	= Potential*eTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double BMagNorm	= BMag*Tau*echarge/MASS;
	double TEMPnorm	= TEMP*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double eDensNorm= eDensity*pow(Radius,3);
	double eTempNorm= eTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double DriftNorm= DriftVel*Tau/(Radius);

	// ***** DEFINE FIELD PARAMETERS ***** //
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.
	threevector BField = BMag*Bhat;
	double DebyeLength 	= sqrt((e0norm*eTempNorm)/eDensNorm);	// Debye Length 
//	double Charge 		= SPEC_CHARGE*PotNorm*(4*PI*e0norm)*exp(-1.0/DebyeLength); 	// Normalised Charge
	double Charge 		= SPEC_CHARGE*PotNorm*(4*PI*e0norm); 	// Normalised Charge,
	double ThermalVel	= sqrt(TEMPnorm);			// Normalised Temperature
	double RhoTherm; 						// Thermal GyroRadius normalised to dust grain radii
        if( BMag <= 0.0 ){
                RhoTherm        = 0.0;
        }else{
                RhoTherm        = ThermalVel/BMagNorm;
        }

//	double VelSIUnits 	= ThermalVel*Radius/(Tau);
//	std::cout << "\nTimeStep = " << TimeStep;
//	std::cout << "\nTEMPnorm = " << TEMPnorm << "\nVelSIUnits = " << VelSIUnits;
//	std::cout << "\nDebyeLength = " << DebyeLength << "\nCharge = " << Charge;
//	std::cout << "\nRhoTherm = " << RhoTherm << "\nRhoThermSI = " << MASS*VelSIUnits/(echarge*BMag*Radius);
//	std::cout << "\nDebyeLength = " << sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2)))/Radius; std::cin.get();


	// ***** DEFINE SIMULATION SPACE ***** //
//	double DebyeHuckelImpactParameter = DebyeLength*LambertW((pow(naturale,2))/DebyeLength);
//	double CoulombImpactParameter	= pow(pow(Charge,2)/(pow(4*PI,2)*e0norm*pow(ThermalVel,2)),0.25);
//	double CoulombImpactParameter	= pow(pow(Charge/(4.0*PI*e0norm),2)/(pow(ThermalVel,2)+BMagNorm/u0norm),0.25);
	double CoulombImpactParameter	= fabs(Charge/(2*PI*e0norm*pow(ThermalVel,2))); // Balance Coulomb to kinetic energy
//	double ImpactParameter = (1.0+2.0*RhoPerp+CoulombImpactParameter);
        double ImpactParameter;
        if( BMag <= 0.0 ){
                ImpactParameter = 1.0+zMaxDebye*CoulombImpactParameter;
        }else{
                ImpactParameter = 1.0+ImpactPar*RhoTherm;
        }
	double zmax 	= 1.0+zMaxDebye*CoulombImpactParameter;	// Top of Simulation Domain, in Dust Radii
	double zmin 	= -1.0-zMinDebye*CoulombImpactParameter;	// Bottom of Simulation Domain, in Dust Radii
//	double zmax 	= 1.0+zMaxDebye*DebyeLength;	// Top of Simulation Domain, in Dust Radii
//	double zmin 	= -1.0-zMinDebye*DebyeLength;	// Bottom of Simulation Domain, in Dust Radii
//	std::cout << "\nDHIP = " << DebyeHuckelImpactParameter;
//	std::cout << "\nRhoTherm = " << RhoTherm;
//	std::cout << "\nCIP = " << CoulombImpactParameter;
//	std::cout << "\nIP = " << ImpactParameter; std::cin.get();


	// ***** DEFINE RANDOM NUMBER GENERATOR ***** //
	std::random_device rd;		// Create Random Device
	std::mt19937 mt(rd());		// Get Random Method
	std::normal_distribution<double> Gaussdist(DriftNorm,sqrt(ThermalVel));
	std::uniform_real_distribution<double> dist(-ImpactParameter, ImpactParameter); // IONS

	// ***** OPEN DATA FILE WITH HEADER ***** //
	time_t now = time(0);		// Get the time of simulation
	char * dt = ctime(&now);
	AngularMomentumDataFile.open(filename + "_AngMom.txt");	
	AngularMomentumDataFile << "## Angular Momentum Data File ##\n";
	AngularMomentumDataFile << "#Date: " << dt;
	AngularMomentumDataFile << "#Input:\timax\tIP\tzmax\tzmin\telec_temp\telec_dens\tion_temp\tRadius\tDensity\tCharge\t\tBMag\tDebye\n";
	AngularMomentumDataFile << "#\t"<<imax<<"\t"<<ImpactParameter<<"\t"<<zmax<<"\t"<<zmin<<"\t"<<eTemp<<"\t\t"<<eDensity<<"\t\t"<<TEMP<<"\t\t"
					<<Radius<<"\t"<<Density<<"\t"<<Charge<<"\t"<<BMag << "\t" << DebyeLength/Radius << "\n";
	AngularMomentumDataFile << "#Output:\tIonNumber\tCollectedIons\tAngularVelocity\tAngularMomentum\n\n";


	// ***** BEGIN LOOP OVER PARTICLE ORBITS ***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	threevector TotalAngularMom(0.0,0.0,0.0);
	DECLARE_LMSUM();
	DECLARE_AMSUM();
	unsigned int j(0), i(0);
	#pragma omp parallel for private(TimeStep) shared(TotalAngularVel) PRIVATE_FILES()
	for( i=0; i < imax; i ++){
		std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
		// ***** RANDOMISE VELOCITY ***** //
		// If the parallel velocity is < vmin, we lose energy which leads to orbits deviating. So we don't include these
		double zvel = Gaussdist(mt);
//		while( zvel >= 0 ){
//			zvel = Gaussdist(mt);
//		}
		threevector Velocity(Gaussdist(mt),Gaussdist(mt),zvel);	// Start with negative z-velocity

		threevector OldVelocity(0.0,0.0,0.0);
		threevector OldPosition(0.0,0.0,0.0);

		// ***** RANDOMISE POSITION ***** //
//		threevector Position(dist(mt),dist(mt),zmax);
		threevector Position(dist(mt),dist(mt),zmax);
		INITIAL_VEL();		// For energy calculations
		INITIAL_POT();		// For energy calculations


		// ***** DO PARTICLE PATH INTEGRATION ***** //
		OPEN_TRACK(filename + std::to_string(i) + ".txt");

		// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
//		threevector EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength,e0norm);
		threevector EField = CoulombField(Position,Charge,e0norm);
		RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);

//		std::cout << "\n" << Position << "\t" << Velocity;
		UpdateVelocityBoris(MASS,EField,BField,-0.5*TimeStep,Velocity);

		// While the particle is not inside the sphere, not outside the simulation domain
		// And not over a specified number of iterations to catch trapped orbits
		int p(0);
		while( Position.mag3() > 1.0 && Position.getz() > zmin && Position.getz() <= zmax && p < 5e5 ){
//			EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength,e0norm);
			EField = CoulombField(Position,Charge,e0norm);
			OldVelocity = Velocity;	// For Angular Momentum Calculations
			OldPosition = Position; // For Angular Momentum Calculations

			TimeStep=(0.001*Position.mag3()/Velocity.mag3()); 	// Adaptive Time Step
			UpdateVelocityBoris(MASS,EField,BField,TimeStep,Velocity);
			Position+=TimeStep*Velocity;

			RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
			p ++;
		}

		if( Position.mag3() < 1.0 ){ // In this case it was captured!
			threevector FinalVelocity = 0.5*(OldVelocity+Velocity);
			threevector FinalPosition = 0.5*(OldPosition+Position);
			threevector CylindricalRadius(FinalPosition.getx(),FinalPosition.gety(),0.0);
			threevector CylindricalVelocity(FinalVelocity.getx(),FinalVelocity.gety(),0.0);
			double DistanceFromAxis = CylindricalRadius.mag3();
			threevector AngularMom = (CylindricalRadius^FinalVelocity); 
			// Moment of Inertia for Hollow Sphere
//			threevector AngularVel = (9*MASS)/(8*PI*NormDens)*	
//				(CylindricalRadius^(FinalVelocity-TotalAngularVel.mag3()*DistanceFromAxis*FinalVelocity.getunit()));
			// Moment of Inertia for solid sphere
//			threevector AngularVel = (15)/(8*PI*NormDens)*	
//				(CylindricalRadius^(FinalVelocity-TotalAngularVel.mag3()*DistanceFromAxis*FinalVelocity.getunit()));

			// SHOULD IT NOT BE:
			threevector AngularVel = (15)/(8*PI*NormDens)*	
				((FinalPosition^FinalVelocity)-(FinalPosition^(TotalAngularVel^FinalPosition)));
			// Moment of Inertia for Solid Sphere, Just Z-Component
//			threevector AngularVel = (15)/(8*PI*NormDens)*
//				((CylindricalRadius^CylindricalVelocity)-(CylindricalRadius^(TotalAngularVel^CylindricalRadius)));

			#pragma omp critical
			{
				TotalAngularVel += AngularVel;
				TotalAngularMom += AngularMom;
				j ++;
				AngularMomentumDataFile << "\n" << i << "\t" << j << "\t" << TotalAngularVel << "\t" << AngularMom;
			}
		}else{
			FINAL_VEL(); 
			FINAL_POT(); 
			OUTPUT_VEL(i); OUTPUT_VEL("\t"); 
			OUTPUT_VEL(100*(FinalVelMag.mag3()/InitialVelMag.mag3()-1));  OUTPUT_VEL("\t");
			OUTPUT_VEL(100*(((FinalVelMag.square()+FinalPot)/(InitialVelMag.square()+InitialPot))-1));  OUTPUT_VEL("\n");
//			OUTPUT_VEL(FinalVelMag.square()-InitialVelMag.square()+FinalPot-InitialPot);  OUTPUT_VEL("\n");
			
			threevector CylindricalRadius(Position.getx(),Position.gety(),0.0);
			LinearMomentumSum += FinalVelMag;
			AngularMomentumSum += CylindricalRadius^FinalVelMag;
		}

		CLOSE_TRACK();
	}
	AngularMomentumDataFile << "\n\n" << LinearMomentumSum << "\t" << AngularMomentumSum;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

	AngularMomentumDataFile.close();
//	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	return 0;
}
