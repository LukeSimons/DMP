#define CALCULATE_MOM
//#define SELF_CONS_CHARGE

//#define SAVE_TRACKS 
#define SAVE_ANGULAR_VEL
//#define SAVE_LINEAR_MOM

//#define TEST_VELPOSDIST
//#define TEST_FINALPOS
//#define TEST_CHARGING
//#define TEST_ANGVEL
//#define TEST_ENERGY

#include <omp.h>	// For parallelisation

#include <iostream>	// for std::cout
#include <array>	
#include <random>	// for std::normal_distribution<> etc.
#include <fstream>	// for std::ofstream
#include <ctime>	// for clock()
#include <math.h> 	// For fabs()
#include <sstream>	// for std::stringstream

//#include "Function.h"	// For LambertW() Function to calculate Debye-Huckel Screening Length
#include "Constants.h"	// Define Pre-processor Directives and constants
#include "threevector.h"// for threevector class
#include "rand_mwts.h"	// Functions to generate correct 1-way maxwellian flux

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\tShow this help message\n\n"
	<< "\t-r,--radius RADIUS\t\t(m), Specify radius of Dust grain\n\n"
	<< "\t-d,--density DENSITY\t\t(m^-^3), Specify density of Dust grain\n\n"
	<< "\t-p,--potential POTENTIAL\t(float), Specify the potential of Dust grain normalised to electron temperature\n\n"
	<< "\t-m,--magfield MAGFIELD\t\t(T), Specify the magnetic field (z direction)\n\n"
	<< "\t-e,--etemp ETEMP\t\t(eV), Specify the temperature of plasma electrons\n\n"
	<< "\t-n,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron density\n\n"
	<< "\t-t,--itemp ITEMP\t\t\t(eV), Specify the temperature of Ions\n\n"
	<< "\t-c,--ichance ICHANCE\t\t\t(eV), Specify the probability of generating an Ion\n\n"
	<< "\t-u,--zmaxcoeff ZMAXCOEFF\t(float), Specify the upper limit of simulation domain as number of distances\n\n"
	<< "\t-l,--zmincoeff ZMINCOEFF\t(float), Specify the lower limit of simulation domain as number of distances\n\n"
	<< "\t-z,--zboundforce ZBOUNDFORCE\t(float), Force the absolute value of simulation domain upper and lower boundaries\n\n"
	<< "\t-b,--impactpar IMPACTPAR\t(float), Specify the radial limit of simulation domain as number of distances\n\n"
	<< "\t-f,--forceimppar FORCEIMPPAR\t(float), Force the absolute value of simulation radial distance\n\n"
	<< "\t-i,--imax IMAX\t(int), Specify the number of particles to be launched\n\n"
	<< "\t-j,--jmax JMAX\t(int), Specify the number of particles to be collected (Over-rides imax)\n\n"
	<< "\t-v,--driftvel DRIFTVEL\t(m s^-^1), Specify the drift velocity of the plasma\n\n"
	<< "\t-s,--seed SEED\t(float), Specify the seed for the random number generator\n\n"
	<< "\t-o,--output OUTPUT\t(string), Specify the suffix of the output file\n\n"
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

void GenerateOrbit(threevector &Position, threevector &Velocity, const double &ImpactParameter, 
			const double &zmin, const double zmax, const double DriftNorm, const double ThermalVel,
			std::mt19937 &mt){

	// ***** DEFINE RANDOM NUMBER GENERATOR ***** //
	std::normal_distribution<double> Gaussdist(DriftNorm,ThermalVel);
	std::uniform_real_distribution<double> rad(0, 1); // IONS

	// ***** RANDOMISE POSITION CYLINDRICALLY ***** //
	double radial_pos=ImpactParameter*sqrt(rad(mt));
	double theta_pos =2*PI*rad(mt);

	Position.setx(radial_pos*cos(theta_pos));
	Position.sety(radial_pos*sin(theta_pos));

	double zvel = rand_mwts(DriftNorm,ThermalVel,mt);	// Generate z-velocity here to determine starting position 
	Position.setz(zmin);
	if( rad(mt) > 0.5 ){
		Position.setz(zmax);
		zvel=-zvel;
	}

	// ***** RANDOMISE VELOCITY ***** //
	Velocity.setx(Gaussdist(mt));
	Velocity.sety(Gaussdist(mt));
	Velocity.setz(zvel);
}

void GyroCentrePos(threevector &Position, threevector &Velocity, double BMagNorm, const threevector &Bhat,
			threevector &GyroCentre2D, double &RhoPerp, const int SPEC_CHARGE){
        threevector PosXY(Position.getx(),Position.gety(),0.0);
        threevector AccelDir = (Velocity.getunit()^Bhat)*SPEC_CHARGE;
        RhoPerp = sqrt(pow(Velocity.getx(),2)+pow(Velocity.gety(),2))/BMagNorm;
        GyroCentre2D = PosXY + (AccelDir*RhoPerp);
}

// THIS HAS BEEN VALIDATED VISUALLY AND WORKS WELL
/*
static void RegenerateMissingOrbits(threevector &Position, threevector &Velocity, const double DriftNorm, const double ThermalVel,
					const double & zmin, const double & zmax, const double &IP, const double & CIP, 
					double BMagNorm, const threevector &Bhat, const double &Charge, 
					unsigned long long & RegeneratedParticles, const int SPEC_CHARGE,
					const double &SpeciesMass, const double &e0norm, std::mt19937 &mt){
	// Calculate orbit parameters and get GyroCentre position
	threevector GyroCentre2D(0.0,0.0,0.0);
	double RhoPerp = 0.0;
	GyroCentrePos(Position,Velocity,BMagNorm,Bhat,GyroCentre2D,RhoPerp,SPEC_CHARGE);
	
	// Check if orbits will intersect the sphere
	if( Charge <= 0 ){	// Repulsive Potential
		// Orbit won't intersect! re-generate orbit
		while( (fabs(GyroCentre2D.mag3() - RhoPerp) > 1.0) ){
//			TRYING TO IMPROVE RE-INJECTION ALGORITHM FOR REPELLED SPECIES
//			HAVING SOME DIFFICULTY AS THESE KIND OF STATEMENTS SEEM TO HAVE NO EFFECT:
//			&& (fabs(Velocity.getz()) < sqrt(fabs(Charge)/(2*PI*SpeciesMass*e0norm))) ){ 
			GenerateOrbit(Position,Velocity,IP,zmin,zmax,DriftNorm,ThermalVel,mt);
			GyroCentrePos(Position,Velocity,BMagNorm,Bhat,GyroCentre2D,RhoPerp,SPEC_CHARGE);
			RegeneratedParticles ++;
		}
	}else if(Charge > 0){	// Attractive Potential
		// Orbit won't intersect! re-generate orbit
		while( fabs(GyroCentre2D.mag3() - RhoPerp) > (1.0 + CIP) ){ 
                        GenerateOrbit(Position,Velocity,IP,zmin,zmax,DriftNorm,ThermalVel,mt);
			GyroCentrePos(Position,Velocity,BMagNorm,Bhat,GyroCentre2D,RhoPerp,SPEC_CHARGE);
                        RegeneratedParticles ++;
                }
	}
}
*/


/*updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation, p.62*/
static void UpdateVelocityBoris(double MASS, threevector Efield, threevector BField, double dt, threevector &Velocity, 
					const int SPEC_CHARGE){

	/*t vector*/
	threevector t = BField*0.5*SPEC_CHARGE*dt;
	
	/*magnitude of t, squared*/
	double t_mag2 = t.square();

	/*s vector*/
	threevector s = 2.0*t*(1/(1+t_mag2));
	
	/*v minus*/
	threevector v_minus = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt; 
	
	/*v prime*/
	threevector v_prime = v_minus + (v_minus^t);
	
	/*v prime*/
	threevector v_plus = v_minus + (v_prime^s);
	
	/*v n+1/2*/
	Velocity = v_plus + Efield*0.5*(SPEC_CHARGE/MASS)*dt;
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
	
	// ***** TIMER AND FILE DECLERATIONS 		***** //
	clock_t begin = clock();
	std::string filename = "Data/DMP";
	std::string suffix	= ".txt";
	DECLARE_TRACK();
	DECLARE_AVEL();			// Data file for angular momentum information
	DECLARE_LMOM();			// Data file for linear momentum information
	std::ofstream RunDataFile;	// Data file for containing the run information

	// ************************************************** //

	// ***** DEFINE DUST PARAMETERS 		***** //
	double Radius 		= 1e-6;	// m, Radius of dust
	double Density 		= 19600;// kg m^-^3, Tungsten
	double Potential	= -2.5;	// Coulombs, Charge in 
	double BMag 		= 1.0; // Tesla, Magnitude of magnetic field

	// ************************************************** //


	// ***** DEFINE PLASMA PARAMETERS 		***** //
	double iTemp 		= 1.0;	// Ion Temperature, eV
	double eTemp 		= 1.0;	// Electron Temperature, eV
	double eDensity		= 1e18;	//1e14;	// m^(-3), Electron density
	double iDensity		= 1e18;	//1e14;	// m^(-3), Ion density
	double DriftVel 	= 0.0;	// m s^-1, This is the Temperature of the species being considered
	double zMaxCoeff	= 3.0;	// Arb, Number of Interaction distances from 0,0 plane to max of simulation domain
	double zMinCoeff	= 3.0;	// Arb, Number of Interaction distances from 0,0 plane to min of simulation domain
	double ZBoundForce	= 0.0;	// Arb, Number of dust grain radii to vertical edge of simulation domain
	double ImpactPar	= 2.0;	// Arb, Multiplicative factor for the Impact Parameter
	double ForceImpPar	= 0.0;	// Arb, Number of dust grain radii to radial edge of simulation domain
	double iChance		= -0.5;	// Arb, Manually set probability of Generating an ion
	unsigned long long imax	= 100;	// Arb, Maximum number of particles to be launched
	unsigned long long jmax	= 2.0e6;// Arb, Number of particles to be collected


	// ************************************************** //


	// ***** RANDOM NUMBER GENERATOR 		***** //
	double seed		= 0.0;	// Arb, Seed for the random number generator
	
	// ************************************************** //


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
		else if( arg == "--itemp" 	|| arg == "-t" )	InputFunction(argc,argv,i,ss0,iTemp);
		else if( arg == "--ichance" 	|| arg == "-c" )	InputFunction(argc,argv,i,ss0,iChance);
		else if( arg == "--zmaxcoeff" 	|| arg == "-u" )	InputFunction(argc,argv,i,ss0,zMaxCoeff);
		else if( arg == "--zmincoeff" 	|| arg == "-l" )	InputFunction(argc,argv,i,ss0,zMinCoeff);
		else if( arg == "--zboundforce"	|| arg == "-z" )	InputFunction(argc,argv,i,ss0,ZBoundForce);
		else if( arg == "--impactpar"	|| arg == "-b" )	InputFunction(argc,argv,i,ss0,ImpactPar);
		else if( arg == "--forceimppar"	|| arg == "-f" )	InputFunction(argc,argv,i,ss0,ForceImpPar);
		else if( arg == "--imax"	|| arg == "-i" )	InputFunction(argc,argv,i,ss0,imax);
		else if( arg == "--jmax"	|| arg == "-j" )	InputFunction(argc,argv,i,ss0,jmax);
		else if( arg == "--driftvel"	|| arg == "-v" )	InputFunction(argc,argv,i,ss0,DriftVel);
		else if( arg == "--seed"	|| arg == "-s" )	InputFunction(argc,argv,i,ss0,seed);
		else if( arg == "--output"	|| arg == "-o" )	InputFunction(argc,argv,i,ss0,suffix);
                else{
			sources.push_back(argv[i]);
		}
	}

	// If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
	double MASS 		= Mp;		// kg, This is the Mass to which quantities are normalised 
	double MassRatio 	= sqrt(Mp/Me);
	double DustMass 	= (4.0/3.0)*PI*pow(Radius,3)*Density;

	// ************************************************** //


	// ***** NORMALISATION 				***** //
	// Normalise TIME to the Gyro-Radius of an Ion at B=100T
	// Normalise MASS to Ion Mass
	// Normalise DISTANCE to Dust Radius
	// Normalise CHARGE to fundamental charge
	double MAGNETIC(100);
        double Tau = MASS/(echarge*100);

	double e0norm 		= epsilon0*MASS*pow(Radius,3)/(pow(echarge*Tau,2));
	double PotNorm		= Potential*eTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));	// NEEDS CHECKING MAYBE
	double DriftNorm	= DriftVel*Tau/(Radius);

	double eDensNorm	= eDensity*pow(Radius,3);
	double iTempNorm	= iTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double eTempNorm	= eTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double DebyeLength 	= sqrt((e0norm*eTempNorm)/eDensNorm);	// Debye Length, THIS WILL BE WRONG FOR ELECTRONS

	// ************************************************** //


	// ***** DEFINE FIELD PARAMETERS 		***** //
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.
	double Charge 		= PotNorm*(4*PI*e0norm); 		// Normalised Charge,

	// ************************************************** //
	

	// ***** DEFINE SIMULATION SPACE 		***** //
	double iThermalVel	= sqrt(iTempNorm/(2.0*PI));		// Normalised Ion Thermal velocity
	double eThermalVel	= sqrt(eTempNorm/(2.0*PI))*MassRatio;	// Normalised Electron Thermal velocity
	double iRhoTherm 	= 0.0;					// Ion gyro-radii are zero by Default
	double eRhoTherm 	= 0.0;					// Electron gyro-radii are zero by Default
	if( BMag != 0.0 ){						// If Magnetic field is non-zero
		// Calculate thermal GyroRadius for ions and electrons normalised to dust grain radii
		iRhoTherm	= iThermalVel/(BMag/MAGNETIC); 
		eRhoTherm	= eThermalVel/(pow(MassRatio,2)*BMag/MAGNETIC); 
	}

	double iCoulombImpactParameter  = fabs(Charge/(2*PI*e0norm*pow(iThermalVel,2))); // Balance Coulomb to kinetic energy
	double eCoulombImpactParameter  = fabs(Charge*pow(MassRatio,2)/(2*PI*e0norm*pow(eThermalVel,2))); // Balance Coulomb to kinetic energy
//	double iImpactParameter = 1.0+ImpactPar*iRhoTherm+zMaxCoeff*iCoulombImpactParameter;
//	double eImpactParameter = 1.0+ImpactPar*eRhoTherm+zMaxCoeff*eCoulombImpactParameter;
	double iImpactParameter = 1.0+ImpactPar*(iRhoTherm+iCoulombImpactParameter); // SUGGESTED CHANGE
	double eImpactParameter = 1.0+ImpactPar*(eRhoTherm+eCoulombImpactParameter); // SUGGESTED CHANGE
	if( ForceImpPar > 0.0 ){
		iImpactParameter = 1.0+ForceImpPar;
		eImpactParameter = 1.0+ForceImpPar;
	}

	double ezmax = 1.05+zMaxCoeff*eCoulombImpactParameter;
	double ezmin = -1.05-zMinCoeff*eCoulombImpactParameter;
	double izmax = 1.0001+zMaxCoeff*iCoulombImpactParameter;
	double izmin = -1.0001-zMinCoeff*iCoulombImpactParameter;
	if( ZBoundForce > 0.0 ){
		ezmax = 1.0+ZBoundForce;
		ezmin = -1.0-ZBoundForce;
		izmax = ezmax;
		izmin = ezmin;
	}

	// ************************************************** //


	// ***** DEFINE PROBABILITY OF ION GENERATION	***** //
	// Define ratio of flux of electrons to ions
//	double ElecToIonRatio = (eDensity/iDensity)*sqrt(eTemp*Mp/(iTemp*Me))*(pow(1+eImpactParameter,2)/pow(1+iImpactParameter,2));
	double ElecToIonRatio = (eDensity/iDensity)*sqrt(eTemp*Mp/(iTemp*Me))*(pow(eImpactParameter,2)/pow(iImpactParameter,2));
	double ProbabilityOfIon = 1.0/(1.0+ElecToIonRatio);
	if( iChance >= 0.0 && iChance <= 1.0 )
		ProbabilityOfIon = iChance;

	// ************************************************** //


	// ***** SEED RANDOM NUMBER GENERATOR IN THREADS***** //
	std::random_device rd;		// Create Random Device
	std::vector<std::mt19937> randnumbers;
	std::uniform_real_distribution<double> rad(0, 1); // IONS
	for(int p = 0; p < omp_get_max_threads(); p ++){
		randnumbers.push_back(std::mt19937(seed+p));
	}

	// ************************************************** //


	// ***** OPEN DATA FILE WITH HEADER 		***** //
	time_t now = time(0);		// Get the time of simulation
	char * dt = ctime(&now);
	OPEN_AVEL();
	OPEN_LMOM();

	RunDataFile.open(filename + suffix);
	RunDataFile << "## Run Data File ##\n";
	RunDataFile << "#Date: " << dt;
	RunDataFile << "#Input:\t\tValue\n\nimax:\t\t"<<imax<<"\njmax:\t\t"<<jmax<<"\nElecToIonratio:\t"<<ElecToIonRatio<<"\nProbOfIon:\t"<<ProbabilityOfIon<<"\n\nElectron Gyro:\t"<<eRhoTherm<<"\nElectron Temp:\t"<<eTemp<<"\nElec Density:\t"<<eDensity<<"\nElectron IP:\t"<<eImpactParameter<<"\nElectron zmax:\t"<<ezmax<<"\nElectron zmin:\t"<<ezmin<<"\n\nIon Gyro:\t"<<iRhoTherm<<"\nIon Temp:\t"<<iTemp<<"\nIon Density:\t"<<iDensity<<"\nIon IP:\t\t"<<iImpactParameter<<"\nIon zmax:\t"<<izmax<<"\nIon zmin:\t"<<izmin<<"\n\nRadius:\t\t"<<Radius<<"\nDensity:\t"<<Density<<"\nCharge:\t\t"<<Charge<<"\nB Field:\t"<<BMag<<"\nDebyeLength:\t"<<DebyeLength/Radius<<"\nDrift Norm:\t"<<DriftNorm<<"\n\n"<<"RNG Seed:\t"<<seed<<"\n\n";

	// ************************************************** //


	// ***** BEGIN LOOP OVER PARTICLE ORBITS 	***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	threevector TotalAngularMom(0.0,0.0,0.0);
	DECLARE_LMSUM();
	DECLARE_AMSUM();

	unsigned long long j(0), i(0), RegeneratedParticles(0), TrappedParticles(0), MissedParticles(0), TotalNum(0);
	long long CapturedCharge(0), RegeneratedCharge(0), TrappedCharge(0), MissedCharge(0), TotalCharge(0);
	#pragma omp parallel for shared(TotalAngularVel,TotalAngularMom,j) PRIVATE_FILES()
	for( i=0; i < imax; i ++){ 	// Loop over maximum number of particles to generate
		if( j <= jmax ){	// Loop until we reach a certain number of particles jmax

//			std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();

			// ***** DETERMINE IF IT'S AN ELECTRON OR ION ***** //
			double BMagNorm = BMag*pow(MassRatio,2)/MAGNETIC;
			double CoulombImpactParameter=eCoulombImpactParameter;
			double ImpactParameter=eImpactParameter;
			double ThermalVel=eThermalVel;
			double zmax= ezmax;        // Top of Simulation Domain, in Dust Radii
			double zmin= ezmin;        // Top of Simulation Domain, in Dust Radii
 			double TimeStep(0.00005);
			double SpeciesMass = 1.0/pow(MassRatio,2);
			int SPEC_CHARGE=-1;
			if( rad(randnumbers[omp_get_thread_num()]) < ProbabilityOfIon ){ // If this is the case, we need to generate an ion
				BMagNorm = BMag/MAGNETIC;
				CoulombImpactParameter=iCoulombImpactParameter;
				ImpactParameter=iImpactParameter;
				ThermalVel=iThermalVel;
				zmax 	= izmax; 
				zmin 	= izmin ;
				TimeStep = 0.005;
				SpeciesMass = 1.0;
				SPEC_CHARGE=1;
			}		
			threevector BField = BMagNorm*Bhat;

			// ************************************************** //
	

			// ***** GENERATE AN ORBIT ***** //
			threevector Position(0.0,0.0,0.0);
			threevector Velocity(0.0,0.0,0.0);
			GenerateOrbit(Position,Velocity,ImpactParameter,zmin,zmax,DriftNorm,ThermalVel,
					randnumbers[omp_get_thread_num()]);

			// ************************************************** //


			// ***** TESTING AREA  				***** //
			// ***** VELOCITY-POSITION DISTRIBUTION TEST 	***** //
			#pragma omp critical 
			{
				PRINT_VPD(Position); PRINT_VPD("\t");	// For debugging look at initial positions
				PRINT_VPD(Velocity); PRINT_VPD("\t");	// For debugging look at initial velocities
				PRINT_VPD( sqrt(pow(Velocity.getx(),2)+pow(Velocity.gety(),2))*SpeciesMass); 
				PRINT_VPD("\n");	// For debugging look at gyro-radii
			}

			// ***** ENERGY TEST: MEASURE INITIAL ENERGY	***** //
			INITIAL_VEL();				// For energy calculations
			INITIAL_POT();				// For energy calculations

			// ***** RECORD TRACK DATA, DEBUG AND TEST	***** //
			OPEN_TRACK(filename + "_Track_" + std::to_string(i) + ".txt");
			RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);

			// ************************************************** //


			// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
//			Calculate Electric Field
//			threevector EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength,e0norm);
			threevector EField = CoulombField(Position,Charge,e0norm);
			UpdateVelocityBoris(SpeciesMass,EField,BField,-0.5*TimeStep,Velocity,SPEC_CHARGE);	

			// ************************************************** //


			// ***** DO PARTICLE PATH INTEGRATION 		***** //
			threevector OldPosition(0.0,0.0,0.0);
			// While we don't exceed a specified number of iterations to catch trapped orbits AND	
			// while the particle is not inside the sphere and not outside the simulation domain
			unsigned int iter(0);
			while( Position.mag3() > 1.0 && Position.getz() >= zmin && Position.getz() <= zmax && iter < 5e5 ){
//				EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength,e0norm);
				EField = CoulombField(Position,Charge,e0norm);
				OldPosition = Position; // For Angular Momentum Calculations

				UpdateVelocityBoris(SpeciesMass,EField,BField,TimeStep,Velocity,SPEC_CHARGE);
				Position+=TimeStep*Velocity;

				RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
				iter ++;
			}	
			CLOSE_TRACK();

			// ************************************************** //


			// ***** PERFORM MOMENTUM CALCULATIONS 		***** //
			threevector FinalPosition = 0.5*(OldPosition+Position);
			threevector AngularMom = SpeciesMass*(FinalPosition^Velocity); 		
			#pragma omp critical
			{
				if( Position.mag3() <= 1.0 ){ // In this case it was captured!
					double AngVelNorm = 5.0*SpeciesMass*MASS/(2.0*DustMass);
					threevector AngularVel = (AngVelNorm)*
					((FinalPosition^Velocity)-(FinalPosition^(TotalAngularVel^FinalPosition)));
					PRINT_FP(fabs(FinalPosition.mag3()-1)); PRINT_FP("\n");
					TotalAngularVel += AngularVel;
					TotalAngularMom += AngularMom;
					j ++;
					CapturedCharge += SPEC_CHARGE;
					PRINT_CHARGE(j)			PRINT_CHARGE("\t")
					PRINT_CHARGE(Charge) 		PRINT_CHARGE("\n")
					PRINT_AVEL((AngVelNorm)*(FinalPosition^Velocity)); PRINT_AVEL("\t");
					PRINT_AVEL((AngVelNorm)*(FinalPosition^Velocity)*(1.0/Tau)); PRINT_AVEL("\n");
					ADD_CHARGE()
					SAVE_AVEL()
					SAVE_LMOM()
				}else if( iter >= 5e5 ){	// In this case it was trapped!
					TrappedParticles ++;
					TrappedCharge += SPEC_CHARGE;
				}else{ 				// In this case it missed!
					LinearMomentumSum += SpeciesMass*Velocity;	
					AngularMomentumSum += AngularMom;
					MissedParticles ++;
					MissedCharge += SPEC_CHARGE;
				} // END OF if ( Position.mag3() < 1.0 )
				FINAL_POT(); 
				PRINT_ENERGY(i); PRINT_ENERGY("\t"); 
				PRINT_ENERGY(100*(Velocity.square()/InitialVel.square()-1.0));  PRINT_ENERGY("\t");
				PRINT_ENERGY(0.5*SpeciesMass*Velocity.square()+SPEC_CHARGE*FinalPot-
						(0.5*SpeciesMass*InitialVel.square()+SPEC_CHARGE*InitialPot));  
				PRINT_ENERGY("\n");
				TotalNum ++;
				TotalCharge += SPEC_CHARGE;
			}
			// ************************************************** //
		}
	} // END OF PARALLELISED FOR LOOP


	// ***** PRINT ANGULAR MOMENTUM AND CHARGE DATA	***** //
	SAVE_MOM("LinMom\t\t\t\tAngMom\n");
	SAVE_MOM(LinearMomentumSum); SAVE_MOM("\t"); SAVE_MOM(AngularMomentumSum); SAVE_MOM("\n\n");
	SAVE_CHARGE("Charge\n")
	SAVE_CHARGE(Charge) SAVE_CHARGE("\n\n")

	// ************************************************** //


	// ***** PRINT CHARGE AND PATH COUNTERS 	***** //
	RunDataFile << "j\tjCharge\tMissed\tMCharge\tRegen\tRCharge\tTrapped\tTCharge\tGross\tGCharge\n"; 
	RunDataFile << j << "\t" << CapturedCharge << "\t" << MissedParticles << "\t" << MissedCharge << "\t" << RegeneratedParticles << "\t" << RegeneratedCharge << "\t" << TrappedParticles << "\t" << TrappedCharge << "\t" << TotalNum << "\t" << TotalCharge << "\n";

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
	RunDataFile << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";

	// ************************************************** //


	// ***** CLOSE DATA FILES 			***** //
	RunDataFile.close();
	CLOSE_AVEL();
	CLOSE_LMOM();

	// ************************************************** //


	return 0;
}
