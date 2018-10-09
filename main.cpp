#define CALCULATE_MOM
//#define CALCULATE_CLOSEST_APPROACH
#define SELF_CONS_CHARGE

//#define SAVE_TRACKS 
#define SAVE_ANGULAR_VEL
#define SAVE_CHARGING
#define SAVE_LINEAR_MOM
//#define SAVE_STARTPOS
#define SAVE_ENDPOS
#define SAVE_SPECIES
//#define SAVE_APPROACH

//#define SPHERICAL_INJECTION
//#define POINT_INJECTION
//#define NO_SPHERE

#define COULOMB_POTENTIAL
//#define DEBYE_POTENTIAL

//#define TEST_VELPOSDIST
//#define TEST_FINALPOS

//#define TEST_CHARGING
//#define TEST_ANGMOM
//#define TEST_ENERGY

#include <omp.h>	// For parallelisation

#include <chrono>	// for chrono::high_resolution_clock::now().time_since_epoch().count();
#include <iostream>	// for std::cout
#include <array>	
#include <random>	// for std::normal_distribution<> etc.
#include <fstream>	// for std::ofstream
#include <ctime>	// for clock()
#include <math.h> 	// for fabs()
#include <sstream>	// for std::stringstream
#include <assert.h>	// for assert()
#include <algorithm>	// for std::min()

//#include "Function.h"	// for LambertW() Function to calculate Debye-Huckel Screening Length
#include "Constants.h"	// Define Pre-processor Directives and constants
#include "threevector.h"// for threevector class
#include "rand_mwts.h"	// Functions to generate correct 1-way maxwellian flux

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\t\tShow this help message\n\n"
	<< "\t-r,--radius RADIUS\t\t(m), Specify radius of Dust grain\n"
	<< "\t\tRadius(=1e-6m) DEFAULT,\t\tBy Default, simulate sphere of size 1um in radius\n\n"
	<< "\t-s,--spin SPIN\t\t(hz), Specify the initial rotation frequency of Dust grain alligned with magnetic field axis\n"
	<< "\t\tSpin(=0.0hz) DEFAULT,\t\tBy Default, simulate sphere with initially zero z angular momentum\n\n"
	<< "\t-a1,--semix SEMIX\t\t(arb), Specify the semi-axis for x in dust radii\n"
	<< "\t\ta1(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
	<< "\t-a2,--semiy SEMIZ\t\t(arb), Specify the semi-axis for y in dust radii\n"
	<< "\t\ta2(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
	<< "\t-a3,--semiz SEMIY\t\t(arb), Specify the semi-axis for z in dust radii\n"
	<< "\t\ta3(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
	<< "\t-d,--density DENSITY\t\t(kgm^-^3), Specify density of Dust grain\n"
	<< "\t\tDensity(=19600kgm^-^3) DEFAULT,\tBy Default, Tungsten Density\n\n"
	<< "\t-p,--potential POTENTIAL\t(double), Specify the potential of Dust grain normalised to electron temperature\n"
	<< "\t\tPotential(=-2.5eV) DEFAULT,\tBy Default, OML Potential in Ti=Te Hydrogen plasma\n\n"
	<< "\t-m,--magfield MAGFIELD\t\t(T), Specify the magnetic field (z direction)\n"
	<< "\t\tBMag(=1.0T) DEFAULT,\t\tBy Default, magnetic field is 1.0T upwards in vertical z\n\n"
	<< "\t-n,--normalised NORMALISED\t(bool), whether normalisation (following Sonmor & Laframboise) is on or off\n"
	<< "\t\tNormalisedVars(=0) DEFAULT,\tFalse, By Default use Tesla and electron volts\n\n"
	<< "\t-te,--etemp ETEMP\t\t(eV), Specify the temperature of plasma electrons\n"
	<< "\t\teTemp(=1eV) DEFAULT,\n\n"
	<< "\t-ne,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron density\n"
	<< "\t\teDensity(=1e18m^-^3) DEFAULT,\tTokamak density\n\n"
	<< "\t-ti,--itemp ITEMP\t\t(eV), Specify the temperature of Ions\n"
	<< "\t\tiTemp(=1eV) DEFAULT,\n\n"
	<< "\t-ni,--idensity IDENSITY\t\t(m^-^3), Specify the plasma ion density\n"
	<< "\t\tiDensity(=1e18m^-^3) DEFAULT,\tTokamak density\n\n"
	<< "\t-c,--ichance ICHANCE\t\t\t(eV), Specify the probability of generating an Ion\n"
	<< "\t\tiChance(=-0.5) DEFAULT,\tFicticious ion generation probability: i.e Self-consistently generate ions & electrons\n\n"
	<< "\t-u,--zmaxcoeff ZMAXCOEFF\t(double), The upper limit of simulation domain as number of Coulomb Interaction lengths\n"
	<< "\t\tzMaxCoeff(=1.0) DEFAULT,\tNumber of Interaction distances from (0,0,Radius) plane to max of simulation domain\n\n"
	<< "\t-l,--zmincoeff ZMINCOEFF\t(double), The lower limit of simulation domain as number of Coulomb Interaction lengths\n"
	<< "\t\tzMaxCoeff(=1.0) DEFAULT,\tNumber of Interaction distances from (0,0,Radius) plane to min of simulation domain\n\n"
	<< "\t-z,--zboundforce ZBOUNDFORCE\t(double), Force the absolute value of simulation domain upper and lower boundaries\n"
	<< "\t\tZBoundForce(=0.0) DEFAULT,\tBy Default, use -u and -l to determine height of injection plane\n\n"
	<< "\t-b,--impactpar IMPACTPAR\t(double), Specify the radial limit of simulation domain as number of distances\n"
	<< "\t\tImpactPar(=2.0) DEFAULT,\tBy Default, Radial extent of injection is three gyro-radii from centre\n\n"
	<< "\t-f,--forceimppar FORCEIMPPAR\t(double), Force the absolute value of simulation radial distance\n"
	<< "\t\tForceImpPar(=0.0) DEFAULT,\tBy Default, use -b to determine radial extent of injection plane\n\n"
	<< "\t-i,--imax IMAX\t\t\t(int), Specify the number of particles to be launched\n"
	<< "\t\timax(=10000) DEFAULT,\t\tBy Default, inject 10,000 particles\n\n"
	<< "\t-j,--jmax JMAX\t\t\t(int), Specify the number of particles to be collected (not exceeding imax)\n"
	<< "\t\tjmax(=5000) DEFAULT,\t\tBy Default, stop simulation if 5,000 particles are collected\n\n"
	<< "\t-rm,--rmax RMAX\t\t\t(int), Specify the number of reflections in z axis before assuming orbit misses\n"
	<< "\t\trmax(=15) DEFAULT,\t\tBy Default, stop simulation if Particle reflects in z axis 15 times\n\n"
	<< "\t-t,--timestep TIMESTEP\t\t(double), Specify the multiplicative factor for the time step\n"
	<< "\t\tTimeStepFactor(=0.0005) DEFAULT,\n\n"
	<< "\t-no,--number NUMBER\t\t\t(int), Specify the number of particles to be captured before saving\n"
	<< "\t\tNumber(=1000) DEFAULT,\t\tBy Default, store average values after collecting 1000 particles\n\n"
	<< "\t-v,--driftvel DRIFTVEL\t\t(m s^-^1), Specify the drift velocity of the plasma\n"
	<< "\t\tDriftVel(=0.0) DEFAULT,\t\tBy Default, No drift velocity\n\n"
	<< "\t-se,--seed SEED\t\t\t(double), Specify the seed for the random number generator\n"
	<< "\t\tSeed(=1.0) DEFAULT,\t\tBy Default, Seed is 1.0\n\n"
	<< "\t-sa,--Saves SAVES\t\t\t(int), Specify the number of particles before saving a run\n"
	<< "\t\tSaves(=1000) DEFAULT,\tBy Default, Save data to Meta datafile once per 1000.0 particles\n\n"
	<< "\t-o,--output OUTPUT\t\t(string), Specify the suffix of the output file\n"
	<< "\t\tsuffix(='.txt') DEFAULT,\tBy Default, Save data to Data/DiMPl.txt\n\n"
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
	// See https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
	// search for "Each component of the velocity vector has a normal distribution"

	std::uniform_real_distribution<double> rad(0.0, 1.0); // IONS
	#ifdef SPHERICAL_INJECTION
		// ***** RANDOMISE POSITION SPHERICALLY ***** //
		double theta_pos  = 2.0*PI*rad(mt);
		double phi_pos    = asin(2.0*rad(mt)-1.0);
		Position.setx(ImpactParameter*cos(phi_pos)*cos(theta_pos));
		Position.sety(ImpactParameter*cos(phi_pos)*sin(theta_pos));
		Position.setz(ImpactParameter*sin(phi_pos));
	
	
		// ***** RANDOMISE VELOCITY SPHERICALLY ***** //
		// http://www.astrosurf.com/jephem/library/li110spherCart_en.htm
		double invel = rand_mwts(0.0,ThermalVel,mt);	// Generate z-velocity here to determine starting position 

		double theta_vel  = 2.0*PI*rad(mt);
                double phi_vel    = asin(2.0*rad(mt)-1.0);

		Velocity.setx(invel*cos(phi_vel)*cos(theta_vel));
		Velocity.sety(invel*cos(phi_vel)*sin(theta_vel));
		Velocity.setz(invel*sin(phi_vel)+DriftNorm);
		

		if( Position*Velocity > 0.0 ){
			Position = -1.0*Position;
		}
	#elif defined(POINT_INJECTION)
		// ***** INJECT AT SINGLE POINT WITH DRIFT VELOCITY ***** //
		Position.setx(0.0);
		Position.sety(0.0);
		Position.setz(zmax);
		
		Velocity.setx(0.0);
		Velocity.sety(0.0);
		Velocity.setz(DriftNorm);
	#else
		// ***** RANDOMISE POSITION CYLINDRICALLY ***** //
		double radial_pos=ImpactParameter*sqrt(rad(mt));
		double theta_pos =2*PI*rad(mt);
		Position.setx(radial_pos*cos(theta_pos));
		Position.sety(radial_pos*sin(theta_pos));
	
		// ***** RANDOMISE VELOCITY CYLINDRICALLY ***** //
	// 	Following work from:
	//	Generating equally weighted test particles from the one-way flux of a drifting Maxwellian
	//	T Makkonen, M I Airila & T Kurki-Suonio
	//	http://iopscience.iop.org/article/10.1088/0031-8949/90/1/015204/data
		double invel = rand_mwts(DriftNorm,ThermalVel,mt);	// Generate z-velocity here to determine starting position 
		Position.setz(zmin);
		if( rad(mt) > 0.5 ){
			Position.setz(zmax);
			invel=-invel+2.0*DriftNorm;
		}
		std::normal_distribution<double> Gaussdist(0.0,ThermalVel);

		Velocity.setx(Gaussdist(mt));
		Velocity.sety(Gaussdist(mt));
		Velocity.setz(invel);
	#endif
}

/*updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation, p.62*/
static void UpdateVelocityBoris(double MASS, const threevector& Efield, const threevector& BField, double dt, 
					threevector &Velocity, const int SPEC_CHARGE){

	/*magnitude of t, squared*/
	double t_mag2 = (BField*0.5*SPEC_CHARGE*dt).square();

	/*v prime*/
	threevector v_prime = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt 
				+ ((Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt)^BField*0.5*SPEC_CHARGE*dt);
	
	/*v prime*/
	threevector v_plus = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt 
				+ (v_prime^(2.0*BField*0.5*SPEC_CHARGE*dt*(1.0/(1.0+t_mag2))));
	
	/*v n+1/2*/
	Velocity = v_plus + Efield*0.5*(SPEC_CHARGE/MASS)*dt;
}

#ifdef COULOMB_POTENTIAL
threevector CoulombField(const threevector &Position, double Charge, double COULOMB_NORM){
//	threevector Efield = (COULOMB_NORM*Charge/Position.square())*Position.getunit(); // 0.07, 0.1475
	return (COULOMB_NORM*Charge/Position.square())*Position.getunit();
}
#endif

#ifdef DEBYE_POTENTIAL
threevector DebyeHuckelField(const threevector &Position, double Charge, double DebyeLength, double COULOMB_NORM){
	if(Charge==0.0) return threevector(0.0,0.0,0.0);
	return COULOMB_NORM*Charge*(1.0/(Position.mag3()))*exp(-(Position.mag3()-1.0)/DebyeLength)
				*(1.0/Position.mag3()+1.0/DebyeLength)*Position.getunit();

//	threevector Efield = COULOMB_NORM*Charge*(1.0/(Position.square()))*exp(-Position.mag3()/DebyeLength)
//				*Position.getunit();
//	return COULOMB_NORM*Charge*(1.0/(Position.square()))*exp(-Position.mag3()/DebyeLength)
//                                *Position.getunit();
}
#endif

int main(int argc, char* argv[]){
	
	// ***** TIMER AND FILE DECLERATIONS 		***** //
	clock_t begin = clock();
	std::string filename = "Data/DiMPl";
	std::string suffix	= ".txt";
	DECLARE_TRACK();
	DECLARE_AVEL();			// Data file for angular momentum information
	DECLARE_LMOM();			// Data file for linear momentum information
	DECLARE_CHA();			// Data file for charge
	DECLARE_SPOS();			// Data file for starting positions
	DECLARE_EPOS();			// Data file for end positions
	DECLARE_APP();			// Date file for closest approaches
	DECLARE_SPEC();			// Data file for the species
	std::ofstream RunDataFile;	// Data file for containing the run information
	std::ofstream InputDataFile;	// Data file for containing the input information

	// ************************************************** //


	// ***** DEFINE DUST PARAMETERS 		***** //
	double Radius 		= 1e-6;		// m, Radius of dust
	double Spin		= 0.0;		// hz, Initial rotation rate
	double Density 		= 19600;	// kg m^-^3, Tungsten
	double Potential	= -2.5;		// Normalised Potential, 
	double BMag 		= 1.0; 		// Tesla, Magnitude of magnetic field
	double BMagIn		= BMag;		// (arb), input magnetic field,in normalised units or Tesla
	bool   NormalisedVars	= false;	// Is magnetic field Normalised according to Sonmor & Laframboise?
	double a1		= 1.0;		// Semi-axis for x in dust-grain radii
	double a2		= 1.0;		// Semi-axis for y in dust-grain radii 
	double a3		= 1.0;		// Semi-axis for z in dust-grain radii

	// ************************************************** //


	// ***** DEFINE PLASMA PARAMETERS 		***** //
	double iTemp 		= 1.0;	// Ion Temperature, eV
	double eTemp 		= 1.0;	// Electron Temperature, eV
	double eDensity		= 1e18;	//1e14;	// m^(-3), Electron density
	double iDensity		= 1e18;	//1e14;	// m^(-3), Ion density
	double DriftVel 	= 0.0;	// m s^-1, This is the Temperature of the species being considered
	double zMaxCoeff	= 1.0;	// Arb, Number of Interaction distances from 0,0 plane to max of simulation domain
	double zMinCoeff	= 1.0;	// Arb, Number of Interaction distances from 0,0 plane to min of simulation domain
	double ZBoundForce	= 0.0;	// Arb, Number of dust grain radii to vertical edge of simulation domain
	double ImpactPar	= 2.0;	// Species gyro-radii, Multiplicative factor for the Impact Parameter
	double ForceImpPar	= 0.0;	// Number of dust grain radii to radial edge of simulation domain
	double iChance		= -0.5;	// Manually set probability of Generating an ion, negative will cause self-consistent
	unsigned long long imax	= 10000;// Arb, Maximum number of particles to be launched
	unsigned long long jmax	= 5000; // Arb, Number of particles to be collected
	unsigned long long num	= 1000; // Arb, Number of particles to be collected before saving
	double TimeStepFactor	= 0.0005;// Arb, Multiplicative factor used to determine size of the timestep
	unsigned int Saves(100);	// Arb, Number of particles to be collected before saving in a run
	unsigned int reflectionsmax(15);// Arb, Number of reflections before rejecting particles


	// ************************************************** //


	// ***** RANDOM NUMBER GENERATOR 		***** //
	// Arb, Seed for the random number generator
	double seed		= 1.0; //std::chrono::high_resolution_clock::now().time_since_epoch().count();
	
	// ************************************************** //


	// ***** DETERMINE USER INPUT ***** //
	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 	|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--radius" 	|| arg == "-r" ) 	InputFunction(argc,argv,i,ss0,Radius);
		else if( arg == "--spin" 	|| arg == "-s" ) 	InputFunction(argc,argv,i,ss0,Spin);
		else if( arg == "--semix" 	|| arg == "-a1") 	InputFunction(argc,argv,i,ss0,a1);
		else if( arg == "--semiy" 	|| arg == "-a2") 	InputFunction(argc,argv,i,ss0,a2);
		else if( arg == "--semiz" 	|| arg == "-a3") 	InputFunction(argc,argv,i,ss0,a3);
		else if( arg == "--density" 	|| arg == "-d" )	InputFunction(argc,argv,i,ss0,Density);
		else if( arg == "--potential" 	|| arg == "-p" )	InputFunction(argc,argv,i,ss0,Potential);
		else if( arg == "--magfield" 	|| arg == "-m" )	InputFunction(argc,argv,i,ss0,BMagIn);
		else if( arg == "--normalised" 	|| arg == "-n" )	InputFunction(argc,argv,i,ss0,NormalisedVars);
		else if( arg == "--etemp" 	|| arg == "-te")	InputFunction(argc,argv,i,ss0,eTemp);
		else if( arg == "--edensity" 	|| arg == "-ne")	InputFunction(argc,argv,i,ss0,eDensity);
		else if( arg == "--itemp" 	|| arg == "-ti")	InputFunction(argc,argv,i,ss0,iTemp);
		else if( arg == "--idensity" 	|| arg == "-ni")	InputFunction(argc,argv,i,ss0,iDensity);
		else if( arg == "--ichance" 	|| arg == "-c" )	InputFunction(argc,argv,i,ss0,iChance);
		else if( arg == "--zmaxcoeff" 	|| arg == "-u" )	InputFunction(argc,argv,i,ss0,zMaxCoeff);
		else if( arg == "--zmincoeff" 	|| arg == "-l" )	InputFunction(argc,argv,i,ss0,zMinCoeff);
		else if( arg == "--zboundforce"	|| arg == "-z" )	InputFunction(argc,argv,i,ss0,ZBoundForce);
		else if( arg == "--impactpar"	|| arg == "-b" )	InputFunction(argc,argv,i,ss0,ImpactPar);
		else if( arg == "--forceimppar"	|| arg == "-f" )	InputFunction(argc,argv,i,ss0,ForceImpPar);
		else if( arg == "--imax"	|| arg == "-i" )	InputFunction(argc,argv,i,ss0,imax);
		else if( arg == "--jmax"	|| arg == "-j" )	InputFunction(argc,argv,i,ss0,jmax);
		else if( arg == "--rmax"	|| arg == "-rm" )	InputFunction(argc,argv,i,ss0,reflectionsmax);
		else if( arg == "--time"	|| arg == "-t" )	InputFunction(argc,argv,i,ss0,TimeStepFactor);
		else if( arg == "--number"	|| arg == "-no")	InputFunction(argc,argv,i,ss0,num);
		else if( arg == "--driftvel"	|| arg == "-v" )	InputFunction(argc,argv,i,ss0,DriftVel);
		else if( arg == "--seed"	|| arg == "-se")	InputFunction(argc,argv,i,ss0,seed);
		else if( arg == "--Saves"	|| arg == "-sa")	InputFunction(argc,argv,i,ss0,Saves);
		else if( arg == "--output"	|| arg == "-o" )	InputFunction(argc,argv,i,ss0,suffix);
                else{
			sources.push_back(argv[i]);
		}
	}
	if( Radius < 0.0 )
		std::cerr << "\nError! Probe Radius is negative\nRadius : " << Radius;
	if( Density < 0.0 )
		std::cerr << "\nError! Probe Density is negative\nDensity : " << Density;
	if( eTemp < 0.0 )
		std::cerr << "\nError! Electron Temperature is negative\neTemp : " << eTemp;
	if( eDensity < 0.0 )
		std::cerr << "\nError! Electron Density is negative\neDensity : " << eDensity;
	if( iTemp < 0.0 )
		std::cerr << "\nError! Ion Temperature is negative\niTemp : " << iTemp;
	if( iDensity < 0.0 )
		std::cerr << "\nError! Ion Density is negative\niDensity : " << iDensity;
	if( iChance < 0.0 )
		std::cout << "\nWarning! Chance of generating Ion is negative, self-consistent flux assumed\niChance : " << iChance;
	if( zMinCoeff < 0.0 )
		std::cerr << "\nError! Lower Vertical Boundary Parameter is negative\nzMinCoeff : " << zMinCoeff;
	if( zMaxCoeff < 0.0 )
		std::cerr << "\nError! Upper Vertical Boundary Parameter is negative\nzMaxCoeff : " << zMaxCoeff;
	if( ZBoundForce < 0.0 )
		std::cerr << "\nError! Force Vertical Boundaries Parameter is negative\nZBoundForce : " << ZBoundForce;
	if( ImpactPar < 0.0 )
		std::cerr << "\nError! Impact Parameter is negative\nImpactPar : " << ImpactPar;
	if( ForceImpPar < 0.0 )
		std::cerr << "\nError! Force Impact Parameter is negative\nForceImpPar : " << ForceImpPar;
	if( imax < jmax )
		std::cerr << "\nError! Total particle goal less than captured particle goal\nimax < jmax : " 
			<< imax << " < " << jmax;
	if( jmax < num )
		std::cout << "\nWarning! Save interval less than captured particle goal. No Angular data recorded\nnum < jmax : " 
			<< num << " < " << jmax;
	if( Saves > 100 && imax <= 100 ){
		std::cout << "\nWarning! Saves greater than number of simulated particles. Code won't run!\nSetting Saves = 1";
		Saves = 1;
	}
		
	InputDataFile.open(filename + "_Input" + suffix);
	InputDataFile << "## Input Data File ##\n";
	InputDataFile << "#Input:\nr\ts\ta1\ta2\ta3\td\tp\tm\tn\tte\tne\tti\tni\tc\tu\tl\tz\tb\tf\ti\tj\tt\tno\tv\tse\tsa\to";
	InputDataFile << "\n" << Radius << "\t" << Spin << "\t" << a1 << "\t" << a2 << "\t" << a3 << "\t" << Density 
		<< "\t" << Potential << "\t" << BMagIn << "\t" << NormalisedVars << "\t" << eTemp << "\t" << eDensity 
		<< "\t" << iTemp << "\t" << iDensity << "\t" << iChance << "\t" << zMaxCoeff << "\t" << zMinCoeff 
		<< "\t" << ZBoundForce << "\t" << ImpactPar << "\t" << ForceImpPar << "\t" << imax << "\t" << jmax 
		<< "\t" << TimeStepFactor << "\t" << num << "\t" << DriftVel << "\t" << seed << "\t" << Saves << "\t" << suffix;
	InputDataFile.close();
	#ifdef SAVE_TRACKS
	if( imax >= 1000 ){
		std::cout << "\nWarning! Attempting to store more than 1000 tracks!\ni = " << imax;
		std::cout << "\nContinue?...\n";
		std::cin.get();
	}
	#endif

	// If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
	double MASS 		= Mp;		// kg, This is the Mass to which quantities are normalised 
	a1 = 1.0/(a1*a1);
	a2 = 1.0/(a2*a2);
	a3 = 1.0/(a3*a3);
	double MassRatio 	= sqrt(Mp/Me);
	double DustMass 	= (4.0/3.0)*PI*pow(Radius,3)*Density;
	if( NormalisedVars ){	// If we're using S&L normalised units but iChance is undertermined

		if( iChance == 0.0 ){ // If we are simulating only Electrons
        	        BMag = sqrt(PI/2.0)*BMagIn*sqrt(Me*eTemp/echarge)/Radius;	// BMag normalised to Electrons
		}else{	// If we are simulating only Ions or otherwise Ions and electrons.
			BMag = sqrt(PI/2.0)*BMagIn*sqrt(Mp*iTemp/echarge)/Radius;	// BMag normalised to Ions
		}

	}else{
		BMag = BMagIn;
                Potential = Potential*echarge/(echarge*eTemp);	// Convert from SI Potential to normalised potential
	}

	// ************************************************** //


	// ***** NORMALISATION 				***** //
	// Normalise TIME to the Gyro-Radius of an Ion at B=100T
	// Normalise MASS to Ion Mass
	// Normalise DISTANCE to Dust Radius
	// Normalise CHARGE to fundamental charge
	double MAGNETIC(100);
        double Tau = MASS/(echarge*MAGNETIC);

	double PotentialNorm	= Potential*(eTemp*echarge)*4*PI*epsilon0*Radius/pow(echarge,2.0);	// Normalised Charge,
	double DriftNorm	= DriftVel*Tau/(Radius);
	double DebyeLength 	= sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2.0)))/Radius;
	double A_Coulomb	= Mp/(4.0*PI*epsilon0*MAGNETIC*MAGNETIC*Radius*Radius*Radius);

	// ************************************************** //


	// ***** DEFINE FIELD PARAMETERS 		***** //
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.

	// ************************************************** //
	

	// ***** DEFINE SIMULATION SPACE 		***** //
	// See : https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
	// For reasoning on choice of Velocity vector
	double iThermalVel	= sqrt(echarge*iTemp/Mp)*(Tau/Radius);		// Normalised Ion Thermal velocity
	double eThermalVel	= sqrt(echarge*eTemp/Me)*(Tau/Radius);	// Normalised Electron Thermal velocity
	double iCoulombImpactParameter  = 10.0*fabs(echarge*echarge*PotentialNorm/(2*PI*epsilon0*echarge*iTemp))/Radius; // Balance Coulomb to kinetic energy
	double eCoulombImpactParameter  = 10.0*fabs(echarge*echarge*PotentialNorm/(2*PI*epsilon0*echarge*eTemp))/Radius; // Balance Coulomb to kinetic energy

	double iRhoTherm = 0.0;					// Ion gyro-radii are zero by Default
	double eRhoTherm = 0.0;					// Electron gyro-radii are zero by Default
//	double TimeStepe = TimeStepFactor/fabs(A_Coulomb*MassRatio*MassRatio); // Normalised time step for electrons
//	double TimeStepi = TimeStepFactor/fabs(A_Coulomb);	// Normalised time step for ions
	double TimeStepe = TimeStepFactor;
	double TimeStepi = TimeStepFactor;
	if( BMag == 0.0 && Potential == 0.0 ){
		TimeStepi = 0.01/iThermalVel;
		TimeStepe = 0.01/eThermalVel;
	}
	if( BMag != 0.0 ){						// If Magnetic field is non-zero
		// Calculate thermal GyroRadius for ions and electrons normalised to dust grain radii
		iRhoTherm	= iThermalVel/(BMag/MAGNETIC); 
		eRhoTherm	= eThermalVel/(pow(MassRatio,2)*BMag/MAGNETIC); 
		
		// Balance electrostatic and magnetic forces to determine height of injection plane
		iCoulombImpactParameter   = sqrt(fabs(echarge*PotentialNorm
						/(4.0*PI*epsilon0*sqrt(echarge*iTemp/(Mp*2.0*PI))*BMag/MAGNETIC)))/Radius;
		eCoulombImpactParameter   = sqrt(fabs(echarge*PotentialNorm
						/(4.0*PI*epsilon0*sqrt(echarge*eTemp/(Me*2.0*PI))*BMag/MAGNETIC)))/Radius;


		if( NormalisedVars ){
			iRhoTherm	= 1.0/BMagIn;
			eRhoTherm	= 1.0/BMagIn;
			
			if( iChance != 0.0 ){ // If there is a finite probability of simulating ions
				eRhoTherm       = 1.0/(BMagIn*MassRatio);
			}
		}
		if( 0.05*eRhoTherm/eThermalVel < TimeStepe || Potential == 0.0 ){
			TimeStepe       = 0.01*eRhoTherm/eThermalVel;   // 1% of Gyro-radius size
		}
		if( 0.05*iRhoTherm/iThermalVel < TimeStepi || Potential == 0.0 ){
	                TimeStepi       = 0.01*iRhoTherm/iThermalVel;   // 1% of Gyro-radius size
		}
	}
	
	// Account for non-spherical impact factor, chose larger of two semi axis
	double semiaxisDistortion=1.0;
	if( a1 < a2 ){
		semiaxisDistortion=a1;
	}else{
		semiaxisDistortion=a2;
	}
	double iImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*(iRhoTherm+iCoulombImpactParameter);
	double eImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*(eRhoTherm+eCoulombImpactParameter);
	if( ForceImpPar > 0.0 ){
		iImpactParameter = 1.0+ForceImpPar;
		eImpactParameter = 1.0+ForceImpPar;
	}

	double ezmax = 1.05/sqrt(a3)+zMaxCoeff*eCoulombImpactParameter;
	double ezmin = -1.05/sqrt(a3)-zMinCoeff*eCoulombImpactParameter;
	double izmax = 1.0001/sqrt(a3)+zMaxCoeff*iCoulombImpactParameter;
	double izmin = -1.0001/sqrt(a3)-zMinCoeff*iCoulombImpactParameter;
	if( ZBoundForce > 0.0 ){
		ezmax = 1.0/sqrt(a3)+ZBoundForce;
		ezmin = -1.0/sqrt(a3)-ZBoundForce;
		izmax = ezmax;
		izmin = ezmin;
	}

	// ************************************************** //

	// ***** DEFINE PROBABILITY OF ION GENERATION	***** //
	// Define ratio of flux of electrons to ions
//	double ElecToIonRatio = (eDensity/iDensity)*sqrt(eTemp*Mp/(iTemp*Me))*(pow(1+eImpactParameter,2)/pow(1+iImpactParameter,2));
	assert(fabs(ezmax)==fabs(ezmin)); // The min and max heights must match as long as the ProbabilityOfIon is the same for
	assert(fabs(izmax)==fabs(izmin)); // both the top and bottom surfaces.
	#ifdef SPHERICAL_INJECTION
		#ifdef COULOMB_POTENTIAL
		double BoltzmanneDensity = eDensity
				*exp(PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*eImpactParameter*Radius*echarge*eTemp));
		double BoltzmanniDensity = iDensity
				*exp(-PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*iImpactParameter*Radius*echarge*iTemp));
		#else
		double BoltzmanneDensity = eDensity
				*exp(PotentialNorm*echarge*echarge*exp(-(eImpactParameter-1.0)/DebyeLength)
				/(4.0*PI*epsilon0*eImpactParameter*Radius*echarge*eTemp));
		double BoltzmanniDensity = iDensity
				*exp(-PotentialNorm*echarge*echarge*exp(-(iImpactParameter-1.0)/DebyeLength)
				/(4.0*PI*epsilon0*iImpactParameter*Radius*echarge*iTemp));
		#endif
	#else
		#ifdef COULOMB_POTENTIAL
		double BoltzmanneDensity = eDensity
				*exp(PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
		double BoltzmanniDensity = iDensity
				*exp(-PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
		#else
		double BoltzmanneDensity = eDensity
                                *exp(PotentialNorm*echarge*echarge*exp(-(ezmax-1.0)/DebyeLength)
				/(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
                double BoltzmanniDensity = iDensity
                                *exp(-PotentialNorm*echarge*echarge*exp(-(izmax-1.0)/DebyeLength)
				/(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
		#endif
	#endif
	double ElecToIonRatio = (BoltzmanneDensity/BoltzmanniDensity)*sqrt(eTemp*Mp/(iTemp*Me))*(pow(eImpactParameter,2)
					/pow(iImpactParameter,2));
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
	OPEN_AVEL();	HEAD_AVEL();
	OPEN_LMOM();	HEAD_LMOM();
	OPEN_CHA();	HEAD_CHA();
	OPEN_EPOS();	HEAD_EPOS();
	OPEN_SPOS();	HEAD_SPOS();
	OPEN_APP();	HEAD_APP()
	OPEN_SPEC();	HEAD_SPEC();

	unsigned long long NumberOfIons = ProbabilityOfIon*imax;
	unsigned long long NumberOfElectrons = (1.0-ProbabilityOfIon)*imax;
	
	RunDataFile.open(filename + suffix);
	RunDataFile << "## Run Data File ##\n";
	RunDataFile << "#Date: " << dt;
	RunDataFile << "#Input:\t\tValue\n\nimax (arb #):\t\t"<<imax<<"\njmax (arb #):\t\t"<<jmax<<"\nIon # (arb #):\t\t" << NumberOfIons<<"\nElec # (arb #):\t\t"<<NumberOfElectrons<<"\nElecToIonratio (arb):\t"<<ElecToIonRatio<<"\nProbOfIon (arb):\t"<<ProbabilityOfIon<<"\n\nElectron Gyro (1/Radius):\t"<<eRhoTherm<<"\nElectron Temp (eV):\t\t"<<eTemp<<"\nElec Density (m^-^3):\t\t"<<BoltzmanneDensity<<"\nElectron IP (1/Radius):\t\t"<<eImpactParameter<<"\nElectron zmax (1/Radius):\t"<<ezmax<<"\nElectron zmin (1/Radius):\t"<<ezmin<<"\nElecs Timestep (Tau):\t\t"<<TimeStepe<<"\n\nIon Gyro (1/Radius):\t"<<iRhoTherm<<"\nIon Temp (eV):\t\t"<<iTemp<<"\nIon Density (m^-^3):\t"<<BoltzmanniDensity<<"\nIon IP (1/Radius):\t"<<iImpactParameter<<"\nIon zmax (1/Radius):\t"<<izmax<<"\nIon zmin (1/Radius):\t"<<izmin<<"\nIon Timestep (Tau):\t"<<TimeStepi<<"\n\nRadius (m):\t\t"<<Radius<<"\nSpin (1/Tau):\t\t"<<Spin*Tau<<"\na1 (1/Radius):\t\t"<<(1.0/sqrt(a1))<<"\na2 (1/Radius):\t\t"<<(1.0/sqrt(a2))<<"\na3 (1/Radius):\t\t"<<(1.0/sqrt(a3))<<"\nDensity (kg m^-^3):\t"<<Density<<"\nCharge (1/echarge):\t\t"<<PotentialNorm<<"\nB Field (T or Radius/GyroRad):\t"<<BMag<<"\nDebyeLength (1/Radius):\t\t"<<DebyeLength <<"\nDrift Norm (Radius/Tau):\t"<<DriftNorm<<"\nTime Norm [Tau] (s):\t\t"<<Tau<<"\n\n"<<"RNG Seed (arb):\t\t"<<seed<<"\nOMP_THREADS (arb):\t"<<omp_get_max_threads()<<"\n\n";
	#ifdef DEBYE_POTENTIAL
		RunDataFile << "* DEBYE HUCKEL POTENTIAL *\n\n";
	#endif
	#ifdef COULOMB_POTENTIAL
		RunDataFile << "* COULOMB POTENTIAL *\n\n";
	#endif

	
	// ************************************************** //


	// ***** BEGIN LOOP OVER PARTICLE ORBITS 	***** //
	threevector TotalAngularVel(0.0,0.0,Spin*Tau);
	threevector TotalAngularMom(0.0,0.0,0.0);
	DECLARE_LMSUM();
	DECLARE_AMSUM();
	DECLARE_AMOM();

	unsigned long long j(0), i(0), RegeneratedParticles(0), TrappedParticles(0), MissedParticles(0), TotalNum(0);
	unsigned long long e_simulated(0), i_simulated(0);
	long long CapturedCharge(0), RegeneratedCharge(0), TrappedCharge(0), MissedCharge(0), TotalCharge(0);

	unsigned long long smax = imax / Saves;

	for( unsigned int s=1; s <= smax; s ++){
		unsigned long long IMAX = (s*imax)/smax;
		RunDataFile.close();
		RunDataFile.clear();
		RunDataFile.open(filename+suffix, std::fstream::app);

		#pragma omp parallel shared(TotalAngularVel,TotalAngularMom,j) PRIVATE_FILES()
		{
		threevector Position(0.0,0.0,0.0);
		threevector Velocity(0.0,0.0,0.0);
		threevector BField(0.0,0.0,0.0);
		threevector EField(0.0,0.0,0.0);
		threevector OldPosition(0.0,0.0,0.0);
		threevector AngularVel(0.0,0.0,0.0);
		threevector InitialPos(0.0,0.0,0.0);
		threevector InitialVel(0.0,0.0,0.0);
		threevector FinalPosition(0.0,0.0,0.0);
		threevector AngularMom(0.0,0.0,0.0);
		unsigned int reflections(0);

		double BMagNorm;
		double ImpactParameter;
		double ThermalVel;
		double zmax;   
		double zmin;   
		double TimeStep;
		double SpeciesMass;

		int SPEC_CHARGE=-1;

		bool EdgeCondition = true;
		bool SphereCondition = true;

		#pragma omp for
		for( i=(IMAX-imax/smax); i < IMAX; i ++){ 	// Loop over maximum number of particles to generate
			if( j <= jmax ){	// Loop until we reach a certain number of particles jmax
	
	//			std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
	
				// ***** DETERMINE IF IT'S AN ELECTRON OR ION ***** //
				// MassRatio is squared because of absence of mass in Boris solver
				#pragma omp critical
				{
				if( rad(randnumbers[omp_get_thread_num()]) < ProbabilityOfIon && i_simulated < NumberOfIons ){ 
					// If this is the case, we need to generate an ion
					BMagNorm = BMag/MAGNETIC;
					ImpactParameter=iImpactParameter;
					ThermalVel=iThermalVel;
					zmax 	= izmax; 
					zmin 	= izmin ;
					TimeStep = TimeStepi;
					SpeciesMass = 1.0;
					SPEC_CHARGE=1;
					i_simulated = i_simulated + 1;
					
				}else{ // If this is the case, we need to generate an electron
					if( e_simulated < NumberOfElectrons ){
						BMagNorm = BMag*pow(MassRatio,2)/MAGNETIC;
						ImpactParameter=eImpactParameter;
						ThermalVel=eThermalVel;
						zmax= ezmax;        // Top of Simulation Domain, in Dust Radii
						zmin= ezmin;        // Top of Simulation Domain, in Dust Radii
		 				TimeStep = TimeStepe;
						SpeciesMass = 1.0/pow(MassRatio,2);
						SPEC_CHARGE=-1;
						e_simulated = e_simulated + 1;
					}else if( i_simulated < NumberOfIons ){
						BMagNorm = BMag/MAGNETIC;
						ImpactParameter=iImpactParameter;
						ThermalVel=iThermalVel;
						zmax 	= izmax; 
						zmin 	= izmin ;
						TimeStep = TimeStepi;
						SpeciesMass = 1.0;
						SPEC_CHARGE=1;
						i_simulated = i_simulated + 1;
					}else{
						std::cerr << "\nWARNING! ALL PARTICLES SIMULATED";
					}
				}
				BField = BMagNorm*Bhat;
				// ***** GENERATE AN ORBIT ***** //

				}
				GenerateOrbit(Position,Velocity,ImpactParameter,zmin,zmax,DriftNorm,ThermalVel,
						randnumbers[omp_get_thread_num()]);

	
				InitialPos = Position;
				InitialVel = Velocity;
				// ************************************************** //
	
				// ***** TESTING AREA  				***** //
				// ***** ANGULAR-MOMENTUM TEST 			***** //
				#ifdef TEST_ANGMOM
				#pragma omp critical 
				{
					ADD_I_AMOM(SpeciesMass*(Position^Velocity));		// For Angular Momentum Calculations
					PRINT_AMOM("IMom = "); PRINT_AMOM(INITIAL_AMOM); PRINT_AMOM("\n");
				}
				#endif
				// ***** VELOCITY-POSITION DISTRIBUTION TEST 	***** //
				#ifdef TEST_VELPOSDIST
				#pragma omp critical
				{
					PRINT_VPD(Position); PRINT_VPD("\t");	// For debugging look at initial positions
					PRINT_VPD(Velocity); PRINT_VPD("\t");	// For debugging look at initial velocities
					PRINT_VPD(Velocity*(Position.getunit())); PRINT_VPD("\t");	
					PRINT_VPD( sqrt(pow(Velocity.getx(),2)+pow(Velocity.gety(),2))*SpeciesMass); 
					PRINT_VPD("\n");	// For debugging look at gyro-radii
				}
				#endif
				// ***** ENERGY TEST: MEASURE INITIAL ENERGY	***** //
				INITIAL_VEL();						// For energy calculations
				INITIAL_POT();						// For energy calculations
				#ifdef CALCULATE_CLOSEST_APPROACH
				double C1 = fabs(2.0*echarge*echarge*PotentialNorm/(Mp*4.0*PI*epsilon0));
				double ri = Position.mag3()*Radius;
				double vmag = Velocity.mag3()*Radius/Tau;
				double vperp = (Position.getunit()^Velocity).mag3()*Radius/Tau;
				double Min_r1 = (C1+sqrt(C1*C1+4.0*(vmag*vmag+C1/ri)*ri*ri*vperp*vperp))/(2.0*(vmag*vmag+C1/ri));
				double Min_r2 = (C1-sqrt(C1*C1+4.0*(vmag*vmag+C1/ri)*ri*ri*vperp*vperp))/(2.0*(vmag*vmag+C1/ri));
				double MinPos = sqrt(zmax*zmax+ImpactParameter*ImpactParameter);
				#endif


				// ***** RECORD TRACK DATA, DEBUG AND TEST	***** //
				OPEN_TRACK(filename + "_Track_" + std::to_string(i) + ".txt");
				RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
	
				// ************************************************** //
	
	
				// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
				// Calculate Electric Field
				#ifdef DEBYE_POTENTIAL
					EField = DebyeHuckelField(Position,PotentialNorm,DebyeLength,A_Coulomb);
				#endif
				#ifdef COULOMB_POTENTIAL
					EField = CoulombField(Position,PotentialNorm,A_Coulomb);
				#endif
				UpdateVelocityBoris(SpeciesMass,EField,BField,-0.5*TimeStep,Velocity,SPEC_CHARGE);	
	
				// ************************************************** //
	
	
				// ***** DO PARTICLE PATH INTEGRATION 		***** //
				OldPosition.setx(0.0);
				OldPosition.sety(0.0);
				OldPosition.setz(0.0);
				// While we don't exceed a specified number of reflections to catch trapped orbits AND	
				// while the particle is not inside the sphere and not outside the simulation domain
				reflections=0;
				unsigned int r_max = reflectionsmax;
				if( PotentialNorm*SPEC_CHARGE > 0.0 ){ // In this case, potential is repulsive
					r_max = 1;
				}

				EdgeCondition = (Position.getz() >= zmin && Position.getz() <= zmax);

				#ifdef SPHERICAL_INJECTION
				EdgeCondition = ((Position.getx()*Position.getx()+Position.gety()*Position.gety()
							+Position.getz()*Position.getz()) <= ImpactParameter*ImpactParameter*1.01);
				#endif
				SphereCondition = (sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2
								+Position.getz()*Position.getz()*a3) > 1.0);
				#ifdef NO_SPHERE
					SphereCondition = true;
				#endif

				while( SphereCondition && EdgeCondition && reflections < r_max ){
					#ifdef DEBYE_POTENTIAL
						EField = DebyeHuckelField(Position,PotentialNorm,DebyeLength,A_Coulomb);
	                                #endif
					#ifdef COULOMB_POTENTIAL
						EField = CoulombField(Position,PotentialNorm,A_Coulomb);
	                                #endif
					OldPosition = Position; // For Angular Momentum Calculations
					double PreviousVelocity = Velocity.getz();
					UpdateVelocityBoris(SpeciesMass,EField,BField,TimeStep,Velocity,SPEC_CHARGE);
					

					if( (Velocity.getz()*PreviousVelocity <= 0.0) ){
						reflections ++;
					}

					Position+=TimeStep*Velocity;
					#ifdef CALCULATE_CLOSEST_APPROACH
					if( OldPosition.mag3() > Position.mag3() || Position.mag3() < MinPos ){
						MinPos = Position.mag3();
					}
					#endif

					EdgeCondition = (Position.getz() >= zmin && Position.getz() <= zmax);
					#ifdef SPHERICAL_INJECTION
	        	                        EdgeCondition = ((Position.getx()*Position.getx()+Position.gety()*Position.gety()
                                        		+Position.getz()*Position.getz()) <= ImpactParameter*ImpactParameter*1.01);
	                                #endif
					SphereCondition = (sqrt(Position.getx()*Position.getx()*a1
							+Position.gety()*Position.gety()*a2
                                                        +Position.getz()*Position.getz()*a3) > 1.0);
				
					#ifdef NO_SPHERE
						SphereCondition = true;
					#endif


					RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
					RECORD_TRACK("\t");
					RECORD_TRACK(sqrt(Position.getx()*Position.getx()
							+Position.gety()*Position.gety()+Position.getz()*Position.getz()));
				}	
				CLOSE_TRACK();
				// ************************************************** //
	
	
				// ***** PERFORM MOMENTUM CALCULATIONS 		***** //
				FinalPosition = 0.5*(OldPosition+Position);
				AngularMom = SpeciesMass*(FinalPosition^Velocity); 		
				#pragma omp critical
				{
					if( sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2+Position.getz()*Position.getz()*a3) < 1.0 ){ // In this case it was captured!
						double AngVelNorm = 5.0*SpeciesMass*MASS/(2.0*DustMass);
						AngularVel = (AngVelNorm)*
						((FinalPosition^Velocity)-(FinalPosition^(TotalAngularVel^FinalPosition)));
//						ADD_F_AMOM(AngularVel*(2.0*DustMass/5.0));


						PRINT_FP(fabs(FinalPosition.mag3()-1)); PRINT_FP("\n");
						TotalAngularVel += AngularVel;
						TotalAngularMom += AngularMom;

						j ++;
						CapturedCharge += SPEC_CHARGE;
						PRINT_CHARGE(j)			PRINT_CHARGE("\t")
						PRINT_CHARGE(PotentialNorm) 	PRINT_CHARGE("\t")
						PRINT_CHARGE(SPEC_CHARGE);	PRINT_CHARGE("\n")
//						PRINT_AMOM((AngVelNorm)*(FinalPosition^Velocity)); PRINT_AMOM("\t");
//						PRINT_AMOM((AngVelNorm)*(FinalPosition^Velocity)*(1.0/Tau)); PRINT_AMOM("\n");
						ADD_CHARGE()
						if(j % num == 0){
							SAVE_AVEL()
							SAVE_LMOM()
							SAVE_CHA()
							SAVE_SPOS()
							SAVE_EPOS()
							SAVE_APP()
							SAVE_SPEC()
						}
						// For SELF_CONS_CHARGE, update the probability of generating ions or electrons
						#ifdef SELF_CONS_CHARGE
						//UPDATE_PROB();
						#endif

					}else if( reflections >= r_max ){        // In this case it was trapped!
						TrappedParticles ++; 
						TrappedCharge += SPEC_CHARGE;
                                        }else{				// In this case it missed!
						SAVE_SPOS()
						SAVE_EPOS()
						SAVE_APP()
						LinearMomentumSum += SpeciesMass*Velocity;	
						AngularMomentumSum += AngularMom;
						MissedParticles ++;
						MissedCharge += SPEC_CHARGE;
					} // END OF if ( Position.mag3() < 1.0 )

					ADD_F_AMOM(SpeciesMass*(Position^Velocity));
					PRINT_AMOM("FMom = "); PRINT_AMOM(FINAL_AMOM); PRINT_AMOM("\n");
					FINAL_POT(); 
					PRINT_ENERGY(i); PRINT_ENERGY("\t"); 
					PRINT_ENERGY(100*(Velocity.square()/InitialVel.square()-1.0));  PRINT_ENERGY("\t");
					PRINT_ENERGY(0.5*Mp*SpeciesMass*InitialVel.square()*Radius*Radius/(Tau*Tau));  PRINT_ENERGY("\t");
					PRINT_ENERGY(0.5*Mp*SpeciesMass*Velocity.square()*Radius*Radius/(Tau*Tau));  PRINT_ENERGY("\t");
					PRINT_ENERGY(SPEC_CHARGE*FinalPot);  PRINT_ENERGY("\t");
					PRINT_ENERGY(SPEC_CHARGE*InitialPot);  PRINT_ENERGY("\t");
					PRINT_ENERGY(0.5*Mp*SpeciesMass*InitialVel.square()*Radius*Radius/(Tau*Tau)
							+SPEC_CHARGE*InitialPot);  PRINT_ENERGY("\t");
					PRINT_ENERGY(0.5*Mp*SpeciesMass*Velocity.square()*Radius*Radius/(Tau*Tau)
							+SPEC_CHARGE*FinalPot);  PRINT_ENERGY("\t");

					PRINT_ENERGY((0.5*Mp*SpeciesMass*Velocity.square()*Radius*Radius/(Tau*Tau)
							+SPEC_CHARGE*FinalPot)/
							(0.5*Mp*SpeciesMass*InitialVel.square()*Radius*Radius/(Tau*Tau)
							+SPEC_CHARGE*InitialPot) - 1.0);  
					PRINT_ENERGY("\n");
					TotalNum ++;
					TotalCharge += SPEC_CHARGE;
				}
			}
		} // END OF PARALLELISED FOR LOOP
		} // END OF PARALLELISED REGION
		// ***** PRINT ANGULAR MOMENTUM AND CHARGE DATA	***** //
		RunDataFile << "*****\n\n\nSave : " << s << "\n\n";
	
		SAVE_MOM("LinMom\t\t\t\tAngMom\n");
		SAVE_MOM(LinearMomentumSum); SAVE_MOM("\t"); SAVE_MOM(AngularMomentumSum); SAVE_MOM("\n\n");
		// Save the Charge in units of electron charges and normalised potential.
		SAVE_CHARGE("Charge\tPotential\n")
		SAVE_CHARGE(PotentialNorm) SAVE_CHARGE("\t") SAVE_CHARGE(PotentialNorm*echarge) 
		SAVE_CHARGE("\n\n")
	
		RunDataFile << "Normalised Ion Current:\tNormalised Electron current\n";
//		Calculate currents for cylindrical geometry with shape factor
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
		/(2.0*(1-cos(asin(sqrt(1.0/(1+izmax*izmax/(iImpactParameter*iImpactParameter))))))*(0.5*(TotalNum+TotalCharge))*iDensity);
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
		/(2.0*(1-cos(asin(sqrt(1.0/(1+ezmax*ezmax/(eImpactParameter*eImpactParameter))))))*(0.5*(TotalNum+TotalCharge))*iDensity) << "\n";
//		Calculate currents for cylindrical geometry
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
					/((TotalNum+TotalCharge)*iDensity);
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
				/((TotalNum-TotalCharge)*eDensity) << "\n";
//		Calculate currents for Spherical geometry with Bfield
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
					/(0.5*(TotalNum+TotalCharge)*iDensity);
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
				/(0.5*(TotalNum-TotalCharge)*eDensity) << "\n";
//		Calculate currents for Spherical geometry without Bfield
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)
					/(0.5*(TotalNum+TotalCharge)*iDensity*(1-cos(asin(1.0/iImpactParameter))));
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)
				/(0.5*(TotalNum-TotalCharge)*eDensity*(1-cos(asin(1.0/eImpactParameter)))) << "\n\n";


		// ************************************************** //
	
	
		// ***** PRINT CHARGE AND PATH COUNTERS 	***** //
		RunDataFile << "j\tjCharge\tMissed\tMCharge\tRegen\tRCharge\tTrapped\tTCharge\tGross\tGCharge\n"; 
		RunDataFile << j << "\t" << CapturedCharge << "\t" << MissedParticles << "\t" << MissedCharge << "\t" << RegeneratedParticles << "\t" << RegeneratedCharge << "\t" << TrappedParticles << "\t" << TrappedCharge << "\t" << TotalNum << "\t" << TotalCharge << "\n";
		if( (i_simulated < NumberOfIons || e_simulated < NumberOfElectrons) && s == smax ){
			std::cerr << "\nError! Total particle goal was not reached! Data may be invalid!";
			RunDataFile << "\n\n* Error! Total particle goal was not reached! *";
			RunDataFile << "\n\n*i_sim =  " << i_simulated << "\te_sim = " << e_simulated << "* ";
		}
		clock_t end = clock();
		double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		RunDataFile << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n\n";
	
		// ************************************************** //
	
	
		// ***** CLOSE DATA FILES 			***** //
		RunDataFile.close();
	}
	CLOSE_AVEL();
	CLOSE_LMOM();
	CLOSE_CHA();
	CLOSE_EPOS();
	CLOSE_SPOS();
	CLOSE_APP();
	CLOSE_SPEC();

	// ************************************************** //


	return 0;
}
