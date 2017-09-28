#ifndef __CONSTANTS_H_INCLUDED__   // if Constants.h hasn't been included yet...
#define __CONSTANTS_H_INCLUDED__


#ifdef STORE_TRACKS 
#define DECLARE_TRACK()	std::ofstream TrackDataFiles	
#define OPEN_TRACK(x)	TrackDataFiles.open(x)
#define RECORD_TRACK(x)	TrackDataFiles << x
#define CLOSE_TRACK()	TrackDataFiles.close()
#else
#define DECLARE_TRACK()
#define OPEN_TRACK(x)
#define RECORD_TRACK(x)
#define CLOSE_TRACK()
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

#endif
