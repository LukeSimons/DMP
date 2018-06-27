#define CALCULATE_MOM
//#define SELF_CONS_CHARGE

//#define SAVE_TRACKS 
#define SAVE_ANGULAR_VEL
#define SAVE_CHARGING
//#define SAVE_LINEAR_MOM
#define SAVE_ENDPOS
#define SAVE_SPECIES

//#define SPHERICAL_INJECTION
//#define POINT_INJECTION
//#define NO_SPHERE

//#define TEST_VELPOSDIST
//#define TEST_FINALPOS
#define TEST_CLOSEST_APPROACH
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

//#include "Function.h"	// for LambertW() Function to calculate Debye-Huckel Screening Length
#include "Constants.h"	// Define Pre-processor Directives and constants
#include "threevector.h"// for threevector class
#include "rand_mwts.h"	// Functions to generate correct 1-way maxwellian flux

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\t\tShow this help message\n\n"
	<< "\t-r,--radius RADIUS\t\t(m), Specify radius of Dust grain\n\n"
	<< "\t-a1,--semix SEMIX\t\t(arb), Specify the semi-axis for x in dust radii\n"
	<< "\t-a2,--semiy SEMIZ\t\t(arb), Specify the semi-axis for y in dust radii\n"
	<< "\t-a3,--semiz SEMIY\t\t(arb), Specify the semi-axis for z in dust radii\n\n"
	<< "\t-d,--density DENSITY\t\t(m^-^3), Specify density of Dust grain\n\n"
	<< "\t-p,--potential POTENTIAL\t(double), Specify the potential of Dust grain normalised to electron temperature\n\n"
	<< "\t-m,--magfield MAGFIELD\t\t(T), Specify the magnetic field (z direction)\n\n"
	<< "\t-n,--normalised NORMALISED\t(bool), whether normalisation (following Sonmor & Laframboise) is on or off\n\n"
	<< "\t-te,--etemp ETEMP\t\t(eV), Specify the temperature of plasma electrons\n\n"
	<< "\t-ne,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron density\n\n"
	<< "\t-ti,--itemp ITEMP\t\t(eV), Specify the temperature of Ions\n\n"
	<< "\t-ni,--idensity IDENSITY\t\t(m^-^3), Specify the plasma ion density\n\n"
	<< "\t-c,--ichance ICHANCE\t\t(eV), Specify the probability of generating an Ion\n\n"
	<< "\t-u,--zmaxcoeff ZMAXCOEFF\t(double), The upper limit of simulation domain as number of Coulomb Interaction lengths\n\n"
	<< "\t-l,--zmincoeff ZMINCOEFF\t(double), The lower limit of simulation domain as number of Coulomb Interaction lengths\n\n"
	<< "\t-z,--zboundforce ZBOUNDFORCE\t(double), Force the absolute value of simulation domain upper and lower boundaries\n\n"
	<< "\t-b,--impactpar IMPACTPAR\t(double), Specify the radial limit of simulation domain as number of distances\n\n"
	<< "\t-f,--forceimppar FORCEIMPPAR\t(double), Force the absolute value of simulation radial distance\n\n"
	<< "\t-i,--imax IMAX\t\t\t(int), Specify the number of particles to be launched\n\n"
	<< "\t-j,--jmax JMAX\t\t\t(int), Specify the number of particles to be collected (not exceeding imax)\n\n"
	<< "\t-no,--number NUMBER\t\t\t(int), Specify the number of particles to be captured before saving\n\n"
	<< "\t-t,--timestep TIMESTEP\t\t(double), Specify the multiplicative factor for the time step\n\n"
	<< "\t-v,--driftvel DRIFTVEL\t\t(m s^-^1), Specify the drift velocity of the plasma\n\n"
	<< "\t-se,--seed SEED\t\t\t(double), Specify the seed for the random number generator\n\n"
	<< "\t-sa,--saves SAVES\t\t\t(int), Specify the number of saves in a run\n\n"
	<< "\t-o,--output OUTPUT\t\t(string), Specify the suffix of the output file\n\n"
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
	std::uniform_real_distribution<double> rad(0.0, 1.0); // IONS

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
		invel=-invel;
	}
	Velocity.setx(Gaussdist(mt));
	Velocity.sety(Gaussdist(mt));
	Velocity.setz(invel);

	#ifdef SPHERICAL_INJECTION
		// ***** RANDOMISE POSITION SPHERICALLY ***** //
		invel = Gaussdist(mt); 
		radial_pos = ImpactParameter;
		theta_pos  = 2.0*PI*rad(mt);
		double phi_pos    = asin(2.0*rad(mt)-1.0);
		Position.setx(radial_pos*cos(phi_pos)*cos(theta_pos));
		Position.sety(radial_pos*cos(phi_pos)*sin(theta_pos));
		Position.setz(radial_pos*sin(phi_pos));
	//	Position.setx(0.0);
	//	Position.sety(0.0);
	//	Position.setz(radial_pos);
	
	
		// ***** RANDOMISE VELOCITY SPHERICALLY ***** //
		// http://www.astrosurf.com/jephem/library/li110spherCart_en.htm
		double theta_vel = Gaussdist(mt)/radial_pos;
		double phi_vel = Gaussdist(mt)/radial_pos;
	
		Velocity.setx(Gaussdist(mt));
		Velocity.sety(Gaussdist(mt));
		Velocity.setz(Gaussdist(mt));
		
		if( Position*Velocity > 0.0 ){
			Position.setx(-1.0*Position.getx());
			Position.sety(-1.0*Position.gety());
			Position.setz(-1.0*Position.getz());
		}
	#endif
	#ifdef POINT_INJECTION
		// ***** INJECT AT SINGLE POINT WITH DRIFT VELOCITY ***** //
		Position.setx(0.0);
		Position.sety(0.0);
		Position.setz(zmax);
		
		Velocity.setx(0.0);
		Velocity.sety(0.0);
		Velocity.setz(DriftNorm);
	#endif
}


/*updates velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation, p.62*/
static void UpdateVelocityBoris(double MASS, threevector Efield, threevector BField, double dt, threevector &Velocity, 
					const int SPEC_CHARGE){

	/*t vector*/
	threevector t = BField*0.5*SPEC_CHARGE*dt;
	
	/*magnitude of t, squared*/
	double t_mag2 = t.square();

	/*s vector*/
	threevector s = 2.0*t*(1.0/(1.0+t_mag2));
	
	/*v minus*/
	threevector v_minus = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt; 
	
	/*v prime*/
	threevector v_prime = v_minus + (v_minus^t);
	
	/*v prime*/
	threevector v_plus = v_minus + (v_prime^s);
	
	/*v n+1/2*/
	Velocity = v_plus + Efield*0.5*(SPEC_CHARGE/MASS)*dt;
}

threevector CoulombField(threevector Position, double Charge, double COULOMB_NORM){
	threevector Efield = (COULOMB_NORM*Charge/Position.square())*Position.getunit(); // 0.07, 0.1475
	return Efield;
}

threevector DebyeHuckelField(threevector Position, double Charge, double DebyeLength, double COULOMB_NORM){
	if(Charge==0.0) return threevector(0.0,0.0,0.0);
//	threevector Efield = COULOMB_NORM*Charge*(1.0/(Position.mag3()))*exp(-Position.mag3()/DebyeLength)
//				*(1.0/Position.mag3()+1.0/DebyeLength)*Position.getunit();

	threevector Efield = COULOMB_NORM*Charge*(1.0/(Position.square()))*exp(-Position.mag3()/DebyeLength)
				*Position.getunit();
	return Efield;
}

int main(int argc, char* argv[]){
	
	// ***** TIMER AND FILE DECLERATIONS 		***** //
	clock_t begin = clock();
	std::string filename = "Data/DiMPl";
	std::string suffix	= ".txt";
	DECLARE_TRACK();
	DECLARE_AVEL();			// Data file for angular momentum information
	DECLARE_LMOM();			// Data file for linear momentum information
	DECLARE_CHA();			// Data file for charge
	DECLARE_EPOS();			// Data file for end positions
	DECLARE_SPEC();			// Data file for the species
	std::ofstream RunDataFile;	// Data file for containing the run information

	// ************************************************** //


	// ***** DEFINE DUST PARAMETERS 		***** //
	double Radius 		= 1e-6;		// m, Radius of dust
	double Density 		= 19600;	// kg m^-^3, Tungsten
	double Potential	= -2.5;		// Normalised Potential, 
	double BMag 		= 1.0; 		// Tesla, Magnitude of magnetic field
	double BMagIn		= BMag;		// (arb), input magnetic field,in normalised units or Tesla
	bool   NormalisedB	= false;	// Is magnetic field Normalised according to Sonmor & Laframboise?
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
	double TimeStepFactor	= 1.0;	// Arb, Multiplicative factor used to determine size of the timestep
	unsigned int saves(2);		// Arb, Number of saves to be performed in a run
	unsigned int reflectionsmax(15);// Arb, Number of reflections before rejecting particles


	// ************************************************** //


	// ***** RANDOM NUMBER GENERATOR 		***** //
	// Arb, Seed for the random number generator
	double seed		= std::chrono::high_resolution_clock::now().time_since_epoch().count();
	
	// ************************************************** //


	// ***** DETERMINE USER INPUT ***** //
	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 	|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--radius" 	|| arg == "-r" ) 	InputFunction(argc,argv,i,ss0,Radius);
		else if( arg == "--semix" 	|| arg == "-a1") 	InputFunction(argc,argv,i,ss0,a1);
		else if( arg == "--semiy" 	|| arg == "-a2") 	InputFunction(argc,argv,i,ss0,a2);
		else if( arg == "--semiz" 	|| arg == "-a3") 	InputFunction(argc,argv,i,ss0,a3);
		else if( arg == "--density" 	|| arg == "-d" )	InputFunction(argc,argv,i,ss0,Density);
		else if( arg == "--potential" 	|| arg == "-p" )	InputFunction(argc,argv,i,ss0,Potential);
		else if( arg == "--magfield" 	|| arg == "-m" )	InputFunction(argc,argv,i,ss0,BMagIn);
		else if( arg == "--normalised" 	|| arg == "-n" )	InputFunction(argc,argv,i,ss0,NormalisedB);
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
		else if( arg == "--time"	|| arg == "-t" )	InputFunction(argc,argv,i,ss0,TimeStepFactor);
		else if( arg == "--number"	|| arg == "-no")	InputFunction(argc,argv,i,ss0,num);
		else if( arg == "--driftvel"	|| arg == "-v" )	InputFunction(argc,argv,i,ss0,DriftVel);
		else if( arg == "--seed"	|| arg == "-se")	InputFunction(argc,argv,i,ss0,seed);
		else if( arg == "--saves"	|| arg == "-sa")	InputFunction(argc,argv,i,ss0,saves);
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
	if( NormalisedB ){	// If we're using S&L normalised units but iChance is undertermined

		if( iChance == 0.0 ){ // If we are simulating only Electrons
        	        BMag = pow(2.0/PI,2.0)*BMagIn*sqrt(PI*Me*eTemp/(2*echarge))/Radius;	// BMag normalised to Electrons
		}else{	// If we are simulating only Ions or otherwise Ions and electrons.
			BMag = pow(2.0/PI,2.0)*BMagIn*sqrt(PI*Mp*iTemp/(2*echarge))/Radius;	// BMag normalised to Ions
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
	double iThermalVel	= sqrt(echarge*iTemp/(Mp*2.0*PI))*(Tau/Radius);		// Normalised Ion Thermal velocity
	double eThermalVel	= sqrt(echarge*eTemp/(Me*2.0*PI))*(Tau/Radius);	// Normalised Electron Thermal velocity
	double iCoulombImpactParameter  = 10.0*fabs(echarge*echarge*PotentialNorm/(2*PI*epsilon0*echarge*iTemp))/Radius; // Balance Coulomb to kinetic energy
	double eCoulombImpactParameter  = 10.0*fabs(echarge*echarge*PotentialNorm/(2*PI*epsilon0*echarge*eTemp))/Radius; // Balance Coulomb to kinetic energy

	double iRhoTherm = 0.0;					// Ion gyro-radii are zero by Default
	double eRhoTherm = 0.0;					// Electron gyro-radii are zero by Default
	double TimeStepe = TimeStepFactor/fabs(A_Coulomb*MassRatio*MassRatio); // Normalised time step for electrons
	double TimeStepi = TimeStepFactor/fabs(A_Coulomb);	// Normalised time step for ions

	if( BMag != 0.0 ){						// If Magnetic field is non-zero
		// Calculate thermal GyroRadius for ions and electrons normalised to dust grain radii
		iRhoTherm	= iThermalVel/(BMag/MAGNETIC); 
		eRhoTherm	= eThermalVel/(pow(MassRatio,2)*BMag/MAGNETIC); 
		
		// Balance electrostatic and magnetic forces to determine height of injection plane
		iCoulombImpactParameter   = sqrt(fabs(echarge*PotentialNorm
						/(4.0*PI*epsilon0*sqrt(echarge*iTemp/(Mp*2.0*PI))*BMag/MAGNETIC)))/Radius;
		eCoulombImpactParameter   = sqrt(fabs(echarge*PotentialNorm
						/(4.0*PI*epsilon0*sqrt(echarge*eTemp/(Me*2.0*PI))*BMag/MAGNETIC)))/Radius;


		if( NormalisedB ){
			iRhoTherm	= 1.0/BMagIn;
			eRhoTherm	= 1.0/BMagIn;
			
			if( iChance != 0.0 ){ // If there is a finite probability of simulating ions
				eRhoTherm       = 1.0/(BMagIn*MassRatio);
			}
		}
		if( 0.05*eRhoTherm/eThermalVel < TimeStepe || Potential == 0.0 ){
			TimeStepe       = 0.05*eRhoTherm/eThermalVel;   // 1% of Gyro-radius size
		}
		if( 0.05*iRhoTherm/iThermalVel < TimeStepi || Potential == 0.0 ){
	                TimeStepi       = 0.05*iRhoTherm/iThermalVel;   // 1% of Gyro-radius size
		}
	}
	
	// Account for non-spherical impact factor, chose larger of two semi axis
	double semiaxisDistortion=1.0;
	if( a1 < a2 ){
		semiaxisDistortion=a1;
	}else{
		semiaxisDistortion=a2;
	}
	double iImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*iRhoTherm+iCoulombImpactParameter;
	double eImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*eRhoTherm+eCoulombImpactParameter;
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
	double BoltzmanneDensity = eDensity*exp(PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
	double BoltzmanniDensity = iDensity*exp(-PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
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
	OPEN_AVEL();
	OPEN_LMOM();
	OPEN_CHA();
	OPEN_EPOS();
	OPEN_SPEC();

	RunDataFile.open(filename + suffix);
	RunDataFile << "## Run Data File ##\n";
	RunDataFile << "#Date: " << dt;
	RunDataFile << "#Input:\t\tValue\n\nimax:\t\t"<<imax<<"\njmax:\t\t"<<jmax<<"\nElecToIonratio:\t"<<ElecToIonRatio<<"\nProbOfIon:\t"<<ProbabilityOfIon<<"\n\nElectron Gyro:\t"<<eRhoTherm<<"\nElectron Temp:\t"<<eTemp<<"\nElec Density:\t"<<BoltzmanneDensity<<"\nElectron IP:\t"<<eImpactParameter<<"\nElectron zmax:\t"<<ezmax<<"\nElectron zmin:\t"<<ezmin<<"\nElecs Timestep:\t"<<TimeStepe<<"\n\nIon Gyro:\t"<<iRhoTherm<<"\nIon Temp:\t"<<iTemp<<"\nIon Density:\t"<<BoltzmanniDensity<<"\nIon IP:\t\t"<<iImpactParameter<<"\nIon zmax:\t"<<izmax<<"\nIon zmin:\t"<<izmin<<"\nIon Timestep:\t"<<TimeStepi<<"\n\nRadius:\t\t"<<Radius<<"\na1:\t\t"<<(1.0/sqrt(a1))<<"\na2:\t\t"<<(1.0/sqrt(a2))<<"\na3:\t\t"<<(1.0/sqrt(a3))<<"\nDensity:\t"<<Density<<"\nCharge:\t\t"<<PotentialNorm<<"\nB Field:\t"<<BMag<<"\nDebyeLength:\t"<<DebyeLength <<"\nDrift Norm:\t"<<DriftNorm<<"\n\n"<<"RNG Seed:\t"<<seed<<"\nOMP_THREADS:\t"<<omp_get_max_threads()<<"\n\n";

	// ************************************************** //


	// ***** BEGIN LOOP OVER PARTICLE ORBITS 	***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	threevector TotalAngularMom(0.0,0.0,0.0);
	DECLARE_LMSUM();
	DECLARE_AMSUM();
	DECLARE_AMOM();

	unsigned long long j(0), i(0), RegeneratedParticles(0), TrappedParticles(0), MissedParticles(0), TotalNum(0);
	long long CapturedCharge(0), RegeneratedCharge(0), TrappedCharge(0), MissedCharge(0), TotalCharge(0);

	for( unsigned int s=1; s <= saves; s ++){
		unsigned long long IMAX = (s*imax)/saves;
		RunDataFile.close();
		RunDataFile.clear();
		RunDataFile.open(filename+suffix, std::fstream::app);
		#pragma omp parallel for shared(TotalAngularVel,TotalAngularMom,j) PRIVATE_FILES()
		for( i=(IMAX-imax/saves); i < IMAX; i ++){ 	// Loop over maximum number of particles to generate
			if( j <= jmax ){	// Loop until we reach a certain number of particles jmax
	
	//			std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
	
				// ***** DETERMINE IF IT'S AN ELECTRON OR ION ***** //
				// MassRatio is squared because of absence of mass in Boris solver
				double BMagNorm = BMag*pow(MassRatio,2)/MAGNETIC;
				double ImpactParameter=eImpactParameter;
				double ThermalVel=eThermalVel;
				double zmax= ezmax;        // Top of Simulation Domain, in Dust Radii
				double zmin= ezmin;        // Top of Simulation Domain, in Dust Radii
	 			double TimeStep(TimeStepe);
				double SpeciesMass = 1.0/pow(MassRatio,2);
				int SPEC_CHARGE=-1;
				if( rad(randnumbers[omp_get_thread_num()]) < ProbabilityOfIon ){ // If this is the case, we need to generate an ion
					BMagNorm = BMag/MAGNETIC;
					ImpactParameter=iImpactParameter;
					ThermalVel=iThermalVel;
					zmax 	= izmax; 
					zmin 	= izmin ;
					TimeStep = TimeStepi;
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
				#ifdef TEST_CLOSEST_APPROACH
				double C1 = fabs(2.0*echarge*echarge*PotentialNorm/(Mp*4.0*PI*epsilon0));
				double ri = Position.mag3()*1e-6;
				double vmag = Velocity.mag3()*1e-6/Tau;
				double vperp = (Position.getunit()^Velocity).mag3()*1e-6/Tau;
				double Min_r1 = (-C1+sqrt(C1*C1+4.0*(vmag*vmag+C1/ri)*ri*ri*vperp*vperp))/(2.0*(vmag*vmag+C1/ri));
				double Min_r2 = (-C1-sqrt(C1*C1+4.0*(vmag*vmag+C1/ri)*ri*ri*vperp*vperp))/(2.0*(vmag*vmag+C1/ri));
//				std::cout << "\nri = " << ri << "\tC1 = " << C1 << "\tvmag = " << vmag << "\tvperp = " << vperp;
//				std::cout << "\nmr1 = " << Min_r1 << "\tmr2 = " << Min_r2 << "\n";
				#pragma omp critical
				{	
				if( Min_r1 <= Radius ){//|| fabs(Min_r2) <= Radius){
					std::cout << "\n" << i << "\t" << j << "\t" << Min_r1 << "\t" << Min_r2;
				}
				}
				#endif


				// ***** RECORD TRACK DATA, DEBUG AND TEST	***** //
				OPEN_TRACK(filename + "_Track_" + std::to_string(i) + ".txt");
				RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
	
				// ************************************************** //
	
	
				// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
				// Calculate Electric Field
				//threevector EField = DebyeHuckelField(Position,PotentialNorm,DebyeLength,A_Coulomb);
				threevector EField = CoulombField(Position,PotentialNorm,A_Coulomb);
				UpdateVelocityBoris(SpeciesMass,EField,BField,-0.5*TimeStep,Velocity,SPEC_CHARGE);	
	
				// ************************************************** //
	
	
				// ***** DO PARTICLE PATH INTEGRATION 		***** //
				threevector OldPosition(0.0,0.0,0.0);
				// While we don't exceed a specified number of reflections to catch trapped orbits AND	
				// while the particle is not inside the sphere and not outside the simulation domain
				unsigned int reflections(0);

				bool EdgeCondition = (Position.getz() >= zmin && Position.getz() <= zmax);

				#ifdef SPHERICAL_INJECTION
				EdgeCondition = ((Position.getx()*Position.getx()+Position.gety()*Position.gety()
							+Position.getz()*Position.getz()) <= ImpactParameter*ImpactParameter*1.01);
				#endif
				bool SphereCondition = (sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2
								+Position.getz()*Position.getz()*a3) > 1.0);
				#ifdef NO_SPHERE
					SphereCondition = true;
				#endif

				while( SphereCondition && EdgeCondition && reflections < reflectionsmax ){
					//EField = DebyeHuckelField(Position,PotentialNorm,DebyeLength,A_Coulomb);
					EField = CoulombField(Position,PotentialNorm,A_Coulomb);
					OldPosition = Position; // For Angular Momentum Calculations

					TimeStep = TimeStepFactor/Velocity.mag3();
					UpdateVelocityBoris(SpeciesMass,EField,BField,TimeStep,Velocity,SPEC_CHARGE);
					
					double PreviousVelocity = Velocity.getz();
					if( (Velocity.getz() > 0.0 && PreviousVelocity <= 0.0) || 
						(Velocity.getz() < 0.0 && PreviousVelocity >= 0.0) ){
						reflections ++;
					}
					Position+=TimeStep*Velocity;

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
				threevector FinalPosition = 0.5*(OldPosition+Position);
				threevector AngularMom = SpeciesMass*(FinalPosition^Velocity); 		
				#pragma omp critical
				{
					if( sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2+Position.getz()*Position.getz()*a3) < 1.0 || reflections >= reflectionsmax ){ // In this case it was captured!
						double AngVelNorm = 5.0*SpeciesMass*MASS/(2.0*DustMass);
						threevector AngularVel = (AngVelNorm)*
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
							SAVE_EPOS()
							SAVE_SPEC()
						}
						// For SELF_CONS_CHARGE, update the probability of generating ions or electrons
						#ifdef SELF_CONS_CHARGE
						UPDATE_PROB();
						#endif
						if( reflections >= reflectionsmax ){        // In this case it was trapped!
        	                                        TrappedParticles ++;
	                                                TrappedCharge += SPEC_CHARGE;
                                        	}
					}else{ 				// In this case it missed!
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
					PRINT_ENERGY(0.5*Mp*SpeciesMass*Velocity.square()+SPEC_CHARGE*FinalPot-
							(0.5*Mp*SpeciesMass*InitialVel.square()+SPEC_CHARGE*InitialPot));  
					PRINT_ENERGY("\n");
					TotalNum ++;
					TotalCharge += SPEC_CHARGE;
				}
			}
		} // END OF PARALLELISED FOR LOOP
		// ***** PRINT ANGULAR MOMENTUM AND CHARGE DATA	***** //
		RunDataFile << "*****\n\n\nSave : " << s << "\n\n";
	
		SAVE_MOM("LinMom\t\t\t\tAngMom\n");
		SAVE_MOM(LinearMomentumSum); SAVE_MOM("\t"); SAVE_MOM(AngularMomentumSum); SAVE_MOM("\n\n");
		// Save the Charge in units of electron charges and normalised potential.
		SAVE_CHARGE("Charge\tPotential\n")
		SAVE_CHARGE(PotentialNorm) SAVE_CHARGE("\t") SAVE_CHARGE(PotentialNorm*echarge) 
		SAVE_CHARGE("\n\n")
	
		RunDataFile << "Normalised Ion Current:\tNormalised Electron current\n";
//		Calculate currents for cylindrical geometry
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
		/(2.0*(1-cos(asin(sqrt(1.0/(1+izmax*izmax/(iImpactParameter*iImpactParameter))))))*(0.5*(TotalNum+TotalCharge))*iDensity);
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
		/(2.0*(1-cos(asin(sqrt(1.0/(1+ezmax*ezmax/(eImpactParameter*eImpactParameter))))))*(0.5*(TotalNum+TotalCharge))*iDensity) << "\n\n";

//		Calculate currents for Spherical geometry with Bfield
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
					/(0.5*(TotalNum+TotalCharge)*iDensity);
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
				/(0.5*(TotalNum-TotalCharge)*eDensity) << "\n\n";
//		Calculate currents for Spherical geometry without Bfield
		RunDataFile << 0.5*BoltzmanniDensity*(j+CapturedCharge)
					/(0.5*(TotalNum+TotalCharge)*iDensity*(1-cos(asin(1.0/iImpactParameter))));
		RunDataFile << "\t\t" << 0.5*BoltzmanneDensity*(j-CapturedCharge)
				/(0.5*(TotalNum-TotalCharge)*eDensity*(1-cos(asin(1.0/eImpactParameter)))) << "\n\n";


		// ************************************************** //
	
	
		// ***** PRINT CHARGE AND PATH COUNTERS 	***** //
		RunDataFile << "j\tjCharge\tMissed\tMCharge\tRegen\tRCharge\tTrapped\tTCharge\tGross\tGCharge\n"; 
		RunDataFile << j << "\t" << CapturedCharge << "\t" << MissedParticles << "\t" << MissedCharge << "\t" << RegeneratedParticles << "\t" << RegeneratedCharge << "\t" << TrappedParticles << "\t" << TrappedCharge << "\t" << TotalNum << "\t" << TotalCharge << "\n";
	
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
	CLOSE_SPEC();

	// ************************************************** //


	return 0;
}
