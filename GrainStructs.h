// Structure containing all the intrinsic and extrinsic information about a material
// Comprised of a particular element

#ifndef __GRAINSTRUCTS_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __GRAINSTRUCTS_H_INCLUDED__

#include "threevector.h"


// NOTE THE FOLLOWING MUST BE WELL DEFINED BY CHILDREN OF MATTER UPON INITIALIZATION
// VapourPressure; 		// N/m^2 or Pa, Pascals
// Emissivity;		// Arb, deviation of EM radiation from black body spectrum
// LinearExpansion; 	// Units of m K^-1
// HeatCapacity;		// kJ/(kg·K)	(Constant Pressure!)

// Data that varies with time and is generally a result of the simulation
struct GrainData{
	// Define if the material is a liquid or a gas
	bool Liquid, Gas;
	bool Breakup;			// True if the dust has undergone liquid breakup

	// ******************************** Do vary with Temperature  ******************************** \\
	// Geometric Data
	double UnheatedRadius;		// metres
	double Mass; 			// Kilogrammes,
	double Radius;			// metres
	double SurfaceArea;		// metres^2
	double Volume;			// metres^3
	double Density;			// Kilogrammes / (metre^3)

	// Thermal Data
	double SuperBoilingTemp;	// Kelvin, Super heated boiling temperature
	double Temperature; 		// Kelvin
	double VapourPressure; 		// N/m^2 or Pa, Pascals
	double Emissivity;		// Arb, deviation of EM radiation from black body spectrum
	double LinearExpansion; 	// Units of m K^-1
	double HeatCapacity;		// kJ/(kg·K)	(Constant Pressure!)

	// Charging Data
	double DeltaSec;		// Secondary Electron Emission Yield (Arb)
	double DeltaTherm;		// Thermionic Electron Emission Yield (Arb)
	double Potential;		// Normalised Grain potential (Arb), Potential = -(e*phi) / (kB * Te)
	bool Positive;			// Defines the sign of the grain charge.

	// Force/Motion Data
	threevector DustPosition;	// Dust position (m)
	threevector DustVelocity;	// Dust velocity (I guess this should be normalised to Cs, but for now: (m/s) )
	double RotationalFrequency;     // Rotational Frequency (s^-1), This is the rate of rotation of the dust.

	// Latent Heat
	double FusionEnergy;		// kJ, Energy put into the material to break bonds of solid
	double VapourEnergy;		// kJ, Energy put into the material to break bond of liquid

};

// Input parameters which are specific to the type of matter being tested
struct ElementConsts{
	char Elem;		// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.

	// Don't vary with Temperature
	double MeltingTemp;	// Kelvin
	double BoilingTemp;	// Kelvin
	double BondEnergy; 	// kJ/mol Amount of energy required to break a single atomic bond
	double WorkFunction;	// eV, Amount of energy required to remove an electron
	double HeatTransAir; 	// W/m^2 K 
	double AtomicMass;	// kg/mol
	double LatentFusion;	// kJ/kg, Energy required to melt 1kg of solid
	double LatentVapour; 	// kJ/kg, Energy required to vapourize 1kg of liquid
	double PlasmaFrequency;	// rad s^-1, plasma frequency of electrons in the material
	double SurfaceTension;	// N/m,
	double RTDensity;	// Kilogrammes / (metre^3), Denisty at room temperature
	double ThermConduct; 	// KiloWatts / (Metre * Kelvin) 

};

#endif
