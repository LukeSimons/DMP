#ifndef __CONSTANTS_H_INCLUDED__   // if Constants.h hasn't been included yet...
#define __CONSTANTS_H_INCLUDED__


#ifdef MODEL_DEBUG
#define Mo_Debug(x) std::cout << x
#else
#define Mo_Debug(x) 
#endif 

#ifdef PLASMAGRID_DEBUG
#define P_Debug(x) std::cout << x
#else
#define P_Debug(x) 
#endif 

#ifdef CHARGING_DEBUG
#define C_Debug(x) std::cout << x
#else
#define C_Debug(x) 
#endif 

#ifdef HEATING_DEBUG
#define H_Debug(x) std::cout << x
#else
#define H_Debug(x) 
#endif 

#ifdef HEATING_DEEP_DEBUG
#define H1_Debug(x) std::cout << x
#else
#define H1_Debug(x) 
#endif 

#ifdef FORCE_DEBUG
#define F_Debug(x) std::cout << x
#else
#define F_Debug(x) 
#endif 

#ifdef FORCE_DEEP_DEBUG
#define F1_Debug(x) std::cout << x
#else
#define F1_Debug(x) 
#endif 

#ifdef MATTER_DEBUG
#define M_Debug(x) std::cout << x
#else
#define M_Debug(x) 
#endif 

#ifdef MATTER_DEEP_DEBUG
#define M2_Debug(x) std::cout << x
#else
#define M2_Debug(x) 
#endif 

#ifdef ELEMENT_DEEP_DEBUG
#define E1_Debug(x) std::cout << x
#else
#define E1_Debug(x) 
#endif 

#ifdef ELEMENT_DEBUG
#define E_Debug(x) std::cout << x
#else
#define E_Debug(x) 
#endif 

#ifdef DTOKSU_DEBUG
#define D_Debug(x) std::cout << x
#else
#define D_Debug(x) 
#endif 

#ifdef DTOKSU_DEEP_DEBUG
#define D1_Debug(x) std::cout << x
#else
#define D1_Debug(x) 
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
