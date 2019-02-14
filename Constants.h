#ifndef __CONSTANTS_H_INCLUDED__   // if Constants.h hasn't been included yet...
#define __CONSTANTS_H_INCLUDED__

#ifdef CALCULATE_MOM
#define DECLARE_MOM()	std::ofstream MomentumDataFile;
#define OPEN_MOM()	MomentumDataFile.open("Data/DiMPl_Momentum.txt");
#define HEAD_MOM()	MomentumDataFile << "#Collect num\tSimulated num\tPx\tPy\tPz\tLx\tLy\tLz";
#define SAVE_MOM()	MomentumDataFile << "\n" << j << "\t" << i << "\t" << LinearMomentumSum << "\t" << AngularMomentumSum;
#define DECLARE_LMSUM()	threevector LinearMomentumSum(0.0,0.0,0.0);
#define DECLARE_AMSUM()	threevector AngularMomentumSum(0.0,0.0,0.0);
#define CLOSE_MOM()	AngularDataFile.close();
#else
#define DECLARE_MOM()
#define OPEN_MOM()
#define HEAD_MOM()
#define SAVE_MOM()
#define DECLARE_LMSUM()
#define DECLARE_AMSUM()
#define CLOSE_MOM()
#endif

#ifdef SELF_CONS_CHARGE
#define ADD_CHARGE()	PotentialNorm += SPEC_CHARGE;
#else
#define ADD_CHARGE()
#endif

#ifdef SAVE_TRACKS 
#define DECLARE_TRACK()	std::ofstream TrackDataFiles	
#define PRIVATE_FILES()	private(TrackDataFiles)
#define OPEN_TRACK(x)	TrackDataFiles.open(x)
#define RECORD_TRACK(x)	TrackDataFiles << x
#define CLOSE_TRACK()	TrackDataFiles.close()
#else
#define DECLARE_TRACK()
#define PRIVATE_FILES()
#define OPEN_TRACK(x)
#define RECORD_TRACK(x)
#define CLOSE_TRACK()
#endif 

#ifdef SAVE_ANGULAR_VEL
#define DECLARE_AVEL()	std::ofstream AngularDataFile;
#define OPEN_AVEL()	AngularDataFile.open("Data/DiMPl_AngVel.txt");
#define HEAD_AVEL()	AngularDataFile << "#Collect num\tSimulated num\tLx\tLy\tLz";
#define SAVE_AVEL()	AngularDataFile << "\n" << j << "\t" << i << "\t" << TotalAngularVel;
#define CLOSE_AVEL()	AngularDataFile.close();
#else
#define DECLARE_AVEL()
#define HEAD_AVEL()
#define OPEN_AVEL()
#define SAVE_AVEL()
#define CLOSE_AVEL()
#endif 

#ifdef SAVE_LINEAR_MOM
#define DECLARE_LMOM()	std::ofstream LinearDataFile;
#define OPEN_LMOM()	LinearDataFile.open("Data/DiMPl_LinMom.txt");
#define HEAD_LMOM()	LinearDataFile << "#Collect num\tSimulated num\tPx_i\tPy_i\tPz_i\tPx_l\tPy_l\tPz_l";
#define SAVE_LMOM()	LinearDataFile << "\n" << j << "\t" << i << "\t" << TotalInjectedMom << "\t" << TotalLostMom;
#define CLOSE_LMOM()	LinearDataFile.close();
#else
#define DECLARE_LMOM()
#define OPEN_LMOM()
#define HEAD_LMOM()
#define SAVE_LMOM()
#define CLOSE_LMOM()
#endif 

#ifdef SAVE_CHARGING
#define DECLARE_CHA()	std::ofstream ChargeDataFile;
#define OPEN_CHA()	ChargeDataFile.open("Data/DiMPl_Charge.txt");
#define HEAD_CHA()	ChargeDataFile << "#Collect num\tSimulated num\tPotential (1/echarge)";
#define SAVE_CHA()	ChargeDataFile << "\n" << j << "\t" << i << "\t" << PotentialNorm;
#define CLOSE_CHA()	ChargeDataFile.close();
#else
#define DECLARE_CHA()
#define OPEN_CHA()
#define HEAD_CHA()
#define SAVE_CHA()
#define CLOSE_CHA()
#endif 

#ifdef SAVE_STARTPOS
#define DECLARE_SPOS()	std::ofstream StartPosDataFile;
#define OPEN_SPOS()	StartPosDataFile.open("Data/DiMPl_StartPos.txt");
#define HEAD_SPOS()	StartPosDataFile << "#Collect num\tSimulated num\tx\ty\tz\tvx\tvy\tvz\tv·r";
#define SAVE_SPOS()	StartPosDataFile << "\n" << j << "\t" << i << "\t" << InitialPos << "\t" << InitialVel << "\t" << InitialVel*(InitialPos.getunit());
#define CLOSE_SPOS()	StartPosDataFile.close();
#else
#define DECLARE_SPOS()
#define OPEN_SPOS()
#define HEAD_SPOS()
#define SAVE_SPOS()
#define CLOSE_SPOS()
#endif 

#ifdef SAVE_ENDPOS
#define DECLARE_EPOS()	std::ofstream EndPosDataFile;
#define OPEN_EPOS()	EndPosDataFile.open("Data/DiMPl_EndPos.txt");
#define HEAD_EPOS()	EndPosDataFile << "#Collect num\tSimulated num\tx\ty\tz\tvx\tvy\tvz\tv·r";
#define SAVE_EPOS()	EndPosDataFile << "\n" << j << "\t" << i << "\t" << Position << "\t" << Velocity << "\t" << Velocity*(Position.getunit());
#define CLOSE_EPOS()	EndPosDataFile.close();
#else
#define DECLARE_EPOS()
#define OPEN_EPOS()
#define HEAD_EPOS()
#define SAVE_EPOS()
#define CLOSE_EPOS()
#endif

#ifdef SAVE_APPROACH
#define DECLARE_APP()	std::ofstream ApproachDataFile;
#define OPEN_APP()	ApproachDataFile.open("Data/DiMPl_Approach.txt");
#define HEAD_APP()	ApproachDataFile << "#Collect num\tSimulated num\trmin";
#define SAVE_APP()	ApproachDataFile << "\n" << j << "\t" << i << "\t" << MinPos;
#define CLOSE_APP()	ApproachDataFile.close();
#else
#define DECLARE_APP()
#define OPEN_APP()
#define HEAD_APP()
#define SAVE_APP()
#define CLOSE_APP()
#endif  

#ifdef SAVE_SPECIES
#define DECLARE_SPEC()	std::ofstream SpeciesDataFile;
#define OPEN_SPEC()	SpeciesDataFile.open("Data/DiMPl_Species.txt");
#define HEAD_SPEC()	SpeciesDataFile << "#Collect num\tSimulated num\tSpecies Charge(echarge)";
#define SAVE_SPEC()	SpeciesDataFile << "\n" << j << "\t" << i << "\t" << SPEC_CHARGE;
#define CLOSE_SPEC()	SpeciesDataFile.close();
#else
#define DECLARE_SPEC()
#define OPEN_SPEC()
#define HEAD_SPEC()
#define SAVE_SPEC()
#define CLOSE_SPEC()
#endif 

#ifdef SAVE_CURRENTS
#define DECLARE_CURR()	std::ofstream CurrentDataFile;
#define OPEN_CURR()	CurrentDataFile.open("Data/DiMPl_Currents.txt");
#define HEAD_CURR()	CurrentDataFile << "#Collect num\tSimulated num\tCylindrical Ii\tCylindrical Ie\tCylindrical Geo Ii\tCylindrical Geo Ie\tSpherical Ii\tSpherical Ie\tSpherical Geo Ii\tSpherical Geo Ie\t";
#define SAVE_CURR()	CurrentDataFile << "\n" << j << "\t" << i << "\t" << CyliCurr << "\t" << CyleCurr << "\t" << CylGeoiCurr << "\t" << CylGeoeCurr << "\t" << SphiCurr << "\t" << SpheCurr << "\t" << SphGeoiCurr << "\t" << SphGeoeCurr;
#define CLOSE_CURR()	CurrentDataFile.close();
#else
#define DECLARE_CURR()
#define OPEN_CURR()
#define HEAD_CURR()
#define SAVE_CURR()
#define CLOSE_CURR()
#endif 

#ifdef SAVE_TOTALS
#define DECLARE_TOT()	std::ofstream TotalDataFile;
#define OPEN_TOT()	TotalDataFile.open("Data/DiMPl_Totals.txt");
#define HEAD_TOT()	TotalDataFile << "#Collect num\tSimulated num\tjCharge\tMissed\tMCharge\tRegen\tRCharge\tTrapped\tTCharge\tGross\tGCharge";
#define SAVE_TOT()	TotalDataFile << "\n" << j << "\t" << i << "\t" << CapturedCharge << "\t" << MissedParticles << "\t" << MissedCharge << "\t" << RegeneratedParticles << "\t" << RegeneratedCharge << "\t" << TrappedParticles << "\t" << TrappedCharge << "\t" << TotalNum << "\t" << TotalCharge;
#define CLOSE_TOT()	TotalDataFile.close();
#else
#define DECLARE_TOT()
#define OPEN_TOT()
#define HEAD_TOT()
#define SAVE_TOT()
#define CLOSE_TOT()
#endif 

#ifdef TEST_VELPOSDIST
#define PRINT_VPD(x)	std::cout << x;
#else
#define PRINT_VPD(x)
#endif

#ifdef TEST_FINALPOS
#define PRINT_FP(x)	std::cout << x;
#else
#define PRINT_FP(x)
#endif

#ifdef TEST_CHARGING
#define PRINT_CHARGE(x)	std::cout << x;
#else
#define PRINT_CHARGE(x)
#endif

#ifdef TEST_ANGMOM
#define DECLARE_AMOM()	threevector INITIAL_AMOM(0.0,0.0,0.0),FINAL_AMOM(0.0,0.0,0.0);
#define ADD_I_AMOM(x)	INITIAL_AMOM = x;
#define ADD_F_AMOM(x)	FINAL_AMOM = x;
#define PRINT_AMOM(x)	std::cout << x;
#else
#define DECLARE_AMOM()
#define ADD_I_AMOM(x)
#define ADD_F_AMOM(x)
#define PRINT_AMOM(x)
#endif

#if defined TEST_COULOMB_ENERGY 
#define C_INITIAL_VEL()	threevector InitialVel = Velocity;
#define C_INITIAL_POT()	double InitialPot = echarge*echarge*PotentialNorm/(4.0*PI*epsilon0*Position.mag3()*Radius);
#define C_FINAL_POT()	double FinalPot = echarge*echarge*PotentialNorm/(4.0*PI*epsilon0*FinalPosition.mag3()*Radius);
#define PRINT_ENERGY(x)	std::cout << x
#else
#define C_INITIAL_VEL()
#define C_INITIAL_POT()
#define C_FINAL_POT()
#define PRINT_ENERGY(x)
#endif

#if defined TEST_DEBYE_ENERGY 
#define D_INITIAL_VEL()	threevector InitialVel = Velocity;
#define D_INITIAL_POT()	double InitialPot = fabs((PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*Radius*Position.mag3()))*exp((Position.mag3()/DebyeLength)-(1.0/DebyeLength)));
#define D_FINAL_POT()	double FinalPot = fabs((PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*Radius*FinalPosition.mag3()))*exp((FinalPosition.mag3()/DebyeLength)-(1.0/DebyeLength)));
#define PRINT_ENERGY(x)	std::cout << x
#else
#define D_INITIAL_VEL()
#define D_INITIAL_POT()
#define D_FINAL_POT()
#define PRINT_ENERGY(x)
#endif



#ifdef PAUSE
#define Pause(); std::cin.get();
#else
#define Pause();
#endif

extern double Kb;  		// (kg m^2 s^-2 K^1) || (J K^-1)
extern double R;		// https://en.wikipedia.org/wiki/Gas_constant 
extern double echarge; 		// C 
extern double Me;		// kg, mass of electron
extern double Mp;		// kg, mass of ion (And Neutrons)
extern double AvNo; 		// mol^-1
extern double PI;
extern double AMU;	 	// kg, Atomic Mass unit
extern double Richardson;	// A/(metres K^2), Richardson constant
extern double c; 		// m/s, Speed of light
extern double h; 		// m^2 kg / s, Planck's Constant
extern double epsilon0; 	// F/m, vacuum permittivity
extern double Sigma;		// Boltsmann constant: W/(m^2 K^4) 
extern double naturale;		// Mathematical constant

#endif
