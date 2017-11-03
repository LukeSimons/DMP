//#define STORE_TRACKS 
#define VEL_POS_DIST
//#define CALCULATE_ENERGY
#define CALCULATE_MOM
#define SAVE_PARTICLE_MOM

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
	<< "\t-j,--jmax JMAX\t(arb), Specify the number of particles to be collected (Over-rides imax)\n\n"
	<< "\t-v,--driftvel DRIFTVEL\t(m s^-^1), Specify the drift velocity of the plasma\n\n"
	<< "\t-o,--output OUTPUT\t(N/A), Specify the suffix of the output file\n\n"
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
			const double &zmin, const double zmax, const double DriftNorm, const double ThermalVel){

	// ***** DEFINE RANDOM NUMBER GENERATOR ***** //
	std::random_device rd;		// Create Random Device
	std::mt19937 mt(rd());		// Get Random Method
	std::normal_distribution<double> Gaussdist(DriftNorm,ThermalVel);
	std::uniform_real_distribution<double> rad(0, 1); // IONS

	// ***** RANDOMISE POSITION CYLINDRICALLY ***** //
	double radial_pos=ImpactParameter*sqrt(rad(mt));
	double theta_pos =2*PI*rad(mt);
	Position.setx(radial_pos*cos(theta_pos));
	Position.sety(radial_pos*sin(theta_pos));
	double zvel = Gaussdist(mt);	// Generate z-velocity here to determine starting position 
	if( zvel >= 0 ){
		Position.setz(zmin);
	}else{
		Position.setz(zmax);
	}

	// ***** RANDOMISE VELOCITY ***** //
	Velocity.setx(Gaussdist(mt));
	Velocity.sety(Gaussdist(mt));
	Velocity.setz(zvel);	

}

void GyroCentrePos(threevector &Position, threevector &Velocity, double BMagNorm, const threevector &Bhat,
			threevector &GyroCentre2D, double &RhoPerp){
	threevector VPerp(Velocity.getx(),Velocity.gety(),0.0);
        threevector PosXY(Position.getx(),Position.gety(),0.0);
        threevector AccelDir = (Velocity.getunit()^Bhat);
        RhoPerp = VPerp.mag3()/BMagNorm;
        GyroCentre2D = PosXY + (AccelDir*RhoPerp);
}

// THIS HAS BEEN VALIDATED VISUALLY AND WORKS WELL
static void RegenerateMissingOrbits(threevector &Position, threevector &Velocity, const double DriftNorm, const double ThermalVel,
					const double & zmin, const double & zmax, const double &IP, const double & CIP, 
					double BMagNorm, const threevector &Bhat, const double &Charge, 
					unsigned long long & RegeneratedParticles){
	// Calculate orbit parameters and get GyroCentre position
	threevector GyroCentre2D(0.0,0.0,0.0);
	double RhoPerp = 0.0;
	GyroCentrePos(Position,Velocity,BMagNorm,Bhat,GyroCentre2D,RhoPerp);
	
	// Check if orbits will intersect the sphere
	if( Charge == 0 ){
		while( fabs(GyroCentre2D.mag3() - RhoPerp) > 1.0 ){ // Orbit won't intersect! re-generate orbit
			GenerateOrbit(Position,Velocity,IP,zmin,zmax,DriftNorm,ThermalVel);
			GyroCentrePos(Position,Velocity,BMagNorm,Bhat,GyroCentre2D,RhoPerp);
			RegeneratedParticles ++;
		}
	}else{	// THIS WON'T WORK FOR ELECTRONS
		while( fabs(GyroCentre2D.mag3() - RhoPerp) > (1.0 + 0.5*CIP)  ){ // Orbit won't intersect! re-generate orbit
                        GenerateOrbit(Position,Velocity,IP,zmin,zmax,DriftNorm,ThermalVel);
			GyroCentrePos(Position,Velocity,BMagNorm,Bhat,GyroCentre2D,RhoPerp);
                        RegeneratedParticles ++;
                }
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
	std::string filename = "Data/DMP_";
	std::string suffix	= "AngMom.txt";
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
	double iDensity		= 1e18;	//1e14;	// m^(-3), Ion density
	double TEMP 		= 1.0;	// eV, This is the Temperature of the species being considered
	double DriftVel 	= 0.0;	// m s^-1, This is the Temperature of the species being considered
	int SPEC_CHARGE		= 1.0;	// arb, This is the charge of the species, should be +1.0 or -1.0 normally 
	double zMaxDebye	= 3.0;	// Arb, Number of debye distances max of simulation is
	double zMinDebye	= 3.0;	// Arb, Number of debye distances max of simulation is
	double ImpactPar	= 2.0;	// Arb, Multiplicative factor for the Impact Parameter
	unsigned long long imax	= 100;	// Arb, Maximum number of particles to be launched
	unsigned long long jmax	= 0.0;// Arb, Number of particles to be collected


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
		else if( arg == "--jmax"	|| arg == "-j" )	InputFunction(argc,argv,i,ss0,jmax);
		else if( arg == "--driftvel"	|| arg == "-v" )	InputFunction(argc,argv,i,ss0,DriftVel);
		else if( arg == "--output"	|| arg == "-o" )	InputFunction(argc,argv,i,ss0,suffix);
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
	double DustMass		= (4.0/3.0)*PI*pow(Radius,3)*Density;

	// ***** NORMALISATION ***** //
	// Normalise TIME to the ratio of the masses and the average time for a ion to hit the dust
	// Normalise MASS to Species Mass
	// Normalise DISTANCE to Dust Radius
	// Normalise CHARGE to fundamental charge
        double Tau;
	double AngVelNorm;
	if( BMag <= 0.0 ){	// THIS PART NEEDS TO BE VALIDATED STILL
		AngVelNorm = sqrt(echarge*TEMP/MASS)*iDensity*pow(Radius,2);
		Tau = sqrt(echarge*TEMP/MASS)*5.0*MASS*iDensity*pow(Radius,2)/(2.0*DustMass);
	}else{
		AngVelNorm = 5.0*MASS/(2.0*DustMass);
		Tau = MASS/(echarge*BMag);
	}
	double e0norm 	= epsilon0*MASS*pow(Radius,3)/(pow(echarge*Tau,2));
	double u0norm 	= (4.0*PI*10e-7*pow(echarge,2))/(MASS*Radius);
 	double TimeStep	= Tau/100000;
	double PotNorm	= Potential*eTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));	// NEEDS CHECKING MAYBE
	double BMagNorm	= BMag*Tau*echarge/MASS;
	double TEMPnorm	= TEMP*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double eDensNorm= eDensity*pow(Radius,3);
	double eTempNorm= eTemp*echarge*pow(Tau,2)/(MASS*pow(Radius,2));
	double DriftNorm= DriftVel*Tau/(Radius);

	// ***** DEFINE FIELD PARAMETERS ***** //
	threevector Bhat(0.0,0.0,1.0);	// Direction of magnetic field, z dir.
	threevector BField = BMagNorm*Bhat;
	double DebyeLength 	= sqrt((e0norm*eTempNorm)/eDensNorm);	// Debye Length 
//	double Charge 		= SPEC_CHARGE*PotNorm*(4*PI*e0norm)*exp(-1.0/DebyeLength); 	// Normalised Charge
	double Charge 		= SPEC_CHARGE*PotNorm*(4*PI*e0norm); 	// Normalised Charge,
	double ThermalVel	= sqrt(2.0*TEMPnorm/PI);			// Normalised Temperature


	// ***** DEFINE SIMULATION SPACE ***** //
//	double DebyeHuckelImpactParameter = DebyeLength*LambertW((pow(naturale,2))/DebyeLength);
//	double CoulombImpactParameter	= pow(pow(Charge,2)/(pow(4*PI,2)*e0norm*pow(ThermalVel,2)),0.25);
//	double CoulombImpactParameter	= pow(pow(Charge/(4.0*PI*e0norm),2)/(pow(ThermalVel,2)+BMagNorm/u0norm),0.25);
	double CoulombImpactParameter	= fabs(Charge/(2*PI*e0norm*pow(ThermalVel,2))); // Balance Coulomb to kinetic energy
        double ImpactParameter;
        if( BMag <= 0.0 ){
                ImpactParameter = 1.0+zMaxDebye*CoulombImpactParameter;
        }else{
		double RhoTherm = ThermalVel/BMagNorm; // Thermal GyroRadius normalised to dust grain radii
                ImpactParameter = 1.0+ImpactPar*RhoTherm+CoulombImpactParameter;
        }

	double zmax 	= 1.5+zMaxDebye*CoulombImpactParameter;	// Top of Simulation Domain, in Dust Radii
	double zmin 	= -1.5-zMinDebye*CoulombImpactParameter;	// Bottom of Simulation Domain, in Dust Radii
//	double zmax 	= 1.0+zMaxDebye*DebyeLength;	// Top of Simulation Domain, in Dust Radii
//	double zmin 	= -1.0-zMinDebye*DebyeLength;	// Bottom of Simulation Domain, in Dust Radii


	// ***** OPEN DATA FILE WITH HEADER ***** //
	time_t now = time(0);		// Get the time of simulation
	char * dt = ctime(&now);
	AngularMomentumDataFile.open(filename + suffix);	
	AngularMomentumDataFile << "## Angular Momentum Data File ##\n";
	AngularMomentumDataFile << "#Date: " << dt;
	AngularMomentumDataFile << "#Input:\timax\tjmax\tIP\tzmax\tzmin\telec_temp\telec_dens\tion_temp\tRadius\tDensity\tCharge\n#Input:\tBMag\tDebye\tDriftNorm\n";
	AngularMomentumDataFile << "#\t"<<imax<<"\t"<<jmax<<"\t"<<ImpactParameter<<"\t"<<zmax<<"\t"<<zmin<<"\t"<<eTemp<<"\t\t"<<eDensity<<"\t\t"<<TEMP<<"\t\t"<<Radius<<"\t"<<Density<<"\t"<<Charge<<"\n#\t"<<BMag << "\t" << DebyeLength/Radius << "\t" << DriftNorm <<"\n";
	AngularMomentumDataFile << "#Output:\tIonNumber\tCollectedIons\tAngularVelocity\tAngularMomentum\n\n";


	// ***** BEGIN LOOP OVER PARTICLE ORBITS ***** //
	threevector TotalAngularVel(0.0,0.0,0.0);
	threevector TotalAngularMom(0.0,0.0,0.0);
	DECLARE_LMSUM();
	DECLARE_AMSUM();
	unsigned long long j(0), i(0), RegeneratedParticles(0), TrappedParticles(0), MissedParticles(0);
	#pragma omp parallel for private(TimeStep) shared(TotalAngularVel) PRIVATE_FILES()
	for( i=0; i < imax; i ++){
	
		if( j != jmax ){
	
//			std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
			threevector Position(0.0,0.0,0.0);
			threevector Velocity(0.0,0.0,0.0);

			GenerateOrbit(Position,Velocity,ImpactParameter,zmin,zmax,DriftNorm,ThermalVel);
		
			// For eliminating orbits that will definitely miss
			if( BMag > 0 )  
				RegenerateMissingOrbits(Position,Velocity,DriftNorm,ThermalVel,zmin,zmax,ImpactParameter,
							CoulombImpactParameter,BMagNorm,Bhat,Charge,RegeneratedParticles); 

			#pragma omp critical 
			{
				PRINT_VP(Position); PRINT_VP("\t");	// For debugging look at initial vel and positions
				PRINT_VP(Velocity); PRINT_VP("\n");	// For debugging look at initial vel and positions
			}


			INITIAL_VEL();				// For energy calculations
			INITIAL_POT();				// For energy calculations


			OPEN_TRACK(filename + "Track_" + std::to_string(i) + ".txt");

			// ***** DO PARTICLE PATH INTEGRATION ***** //
//			Calculate Electric Field
//			threevector EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength,e0norm);
			threevector EField = CoulombField(Position,Charge,e0norm);

			RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);


			// ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
			UpdateVelocityBoris(MASS,EField,BField,-0.5*TimeStep,Velocity);	


			unsigned int iter(0);
			bool Repeat(true);
			threevector OldVelocity(0.0,0.0,0.0), OldPosition(0.0,0.0,0.0);
			// While we don't exceed a specified number of iterations to catch trapped orbits	
			while( Repeat ){	
				// While the particle is not inside the sphere and not outside the simulation domain
				while( Position.mag3() > 1.0 && Position.getz() >= zmin && Position.getz() <= zmax && iter < 5e5 ){
//					EField = DebyeHuckelField(Position,Charge,Radius,eDensity,eTemp,DebyeLength,e0norm);
					EField = CoulombField(Position,Charge,e0norm);
					OldVelocity = Velocity;	// For Angular Momentum Calculations
					OldPosition = Position; // For Angular Momentum Calculations

					TimeStep=0.1;       // Adaptive Time Step
					UpdateVelocityBoris(MASS,EField,BField,TimeStep,Velocity);
					Position+=TimeStep*Velocity;

					RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
					iter ++;
				}
				if( iter < 5e5 ){
					Repeat = false;
				}else{	// Particle is considered trapped
					CLOSE_TRACK();
					GenerateOrbit(Position,Velocity,ImpactParameter,zmin,zmax,DriftNorm,ThermalVel);
					RegenerateMissingOrbits(Position,Velocity,DriftNorm,ThermalVel,zmin,zmax,ImpactParameter,
							CoulombImpactParameter,BMagNorm,Bhat,Charge,RegeneratedParticles);
					INITIAL_VEL();          // For energy calculations
        			        INITIAL_POT();          // For energy calculations				
					OPEN_TRACK(filename + std::to_string(i) + "_retry.txt");
					EField = CoulombField(Position,Charge,e0norm);
		                	RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
					UpdateVelocityBoris(MASS,EField,BField,-0.5*TimeStep,Velocity);
					iter = 0;
					Repeat = true;
					TrappedParticles ++;
				}	
			}

			threevector FinalVelocity = 0.5*(OldVelocity+Velocity);
			threevector FinalPosition = 0.5*(OldPosition+Position);
			threevector CylindricalRadius(FinalPosition.getx(),FinalPosition.gety(),0.0);
			threevector AngularMom = (CylindricalRadius^FinalVelocity); 

			if( Position.mag3() < 1.0 ){ // In this case it was captured!
				threevector AngularVel = (AngVelNorm)*
				((FinalPosition^FinalVelocity)-(FinalPosition^(TotalAngularVel^FinalPosition)));
				#pragma omp critical
				{
					if( j < jmax){
						TotalAngularVel += AngularVel;
						TotalAngularMom += AngularMom;
						j ++;
						SAVE_MOM()
					}
				}
			}else{
				FINAL_VEL(); 
				FINAL_POT(); 
				OUTPUT_VEL(i); OUTPUT_VEL("\t"); 
				OUTPUT_VEL(100*(FinalVelMag.mag3()/InitialVelMag.mag3()-1));  OUTPUT_VEL("\t");
				OUTPUT_VEL(100*(((FinalVelMag.square()+FinalPot)/(InitialVelMag.square()+InitialPot))-1));  
				OUTPUT_VEL("\n");
			
				LinearMomentumSum += FinalVelocity;
				AngularMomentumSum += AngularMom;
				MissedParticles ++;
			} // END OF if ( Position.mag3() < 1.0 )
		} // END OF if (j != jmax-1)
		CLOSE_TRACK();
	} // END OF PARALLELISED FOR LOOP
	AngularMomentumDataFile << "\n\n" << LinearMomentumSum << "\t" << AngularMomentumSum << "\t" 
		<< j << "\t" << MissedParticles << "\t" << RegeneratedParticles << "\t" << TrappedParticles;

	AngularMomentumDataFile.close();

//	clock_t end = clock();
//	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
//	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	return 0;
}
