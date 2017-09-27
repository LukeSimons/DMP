//#define ELEMENT_DEBUG
//#define ELEMENT_DEEP_DEBUG

#include "Tungsten.h"
#include "Constants.h"

const struct ElementConsts TungstenConsts = {
	'W',			// Specifies the element
	3422, 			// K, Melting temperature at atmospheric pressure
	5555,			// K, Boiling temperature at atmospheric pressure
	774,			// Units of kJ/mol, Bond Energy, Estimated as being equal to Latent Vapour Energy
	3.4,			// ev, Work Function, Taken from DTOKS and matched with wikipedia
	5.7,	 		// kJ/(m^2 K), HeatTransfer coefficient, 
					// http://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
	0.18384,		// AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
	35.3/0.18384, 		// kJ/kg, Latent Fusion Energy From Wikipedia
	774/0.18384, 		// kJ/kg, Latent Vapour Energy From Wikipedia
	9.12e15,		// rad s^-1, Plasma Frequency, 'On thermal radiation from heated metallic dust grains' a. Pigarov
	2.333,			// N/m, Surface Tension B. Keene, 1993, Review of data for the surface tension of pure metals
	19600,			// (kg/m^3) from Wikipedia, Room temperature density
	0.163			// kW/m K at 20 degrees celsius
};

Tungsten::Tungsten():Matter(&TungstenConsts){
	E_Debug("\n\nIn Tungsten::Tungsten():Matter(&TungstenConsts)\n\n");
	tungsten_defaults();
	update();
}

Tungsten::Tungsten(double radius):Matter(radius,&TungstenConsts){
	E_Debug("\n\nIn Tungsten::Tungsten(double radius):Matter(radius,&TungstenConsts)\n\n");
	tungsten_defaults();
	update();
}

Tungsten::Tungsten(double radius, double tempin):Matter(radius,tempin,&TungstenConsts){
	E_Debug("\n\nIn Tungsten::Tungsten(double radius, double tempin):Matter(radius,tempin,&TungstenConsts)\n\n");
	tungsten_defaults();
	update_state(0.0);		// Temperature dependant
	update_models('c','c','c','y');
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Tungsten::Tungsten(double radius, double tempin, std::array<char,4> &constmodels):Matter(radius,tempin,&TungstenConsts){
	E_Debug("\n\nIn Tungsten::Tungsten(double radius, double tempin, std::array<char,4> &constmodels):Matter(radius,tempin,&TungstenConsts)\n\n\t");
	tungsten_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Tungsten::tungsten_defaults(){
	E_Debug("\tIn Tungsten::tungsten_defaults()\n\n");

	// http://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.html
        St.HeatCapacity = 0.13398; 		// kJ/(kg-K)         
	St.Emissivity = 0.04; 			// Arb, http://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
        St.SuperBoilingTemp = Ec.BoilingTemp; 	// K, At any pressure
	St.Density = 19600;			// (kg/m^3) from Wikipedia
	update_models('c','c','c','y');
	

//	St.ThermConduct = 0.163;		// kW/m K at 20 degrees celsius
}

/// ************************************************************************************************* \\\

void Tungsten::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\tIn Tungsten::update_heatcapacity()\n\n");

	if( St.Temperature < 300 ){ 

		static bool runOnce = true;
		WarnOnce(runOnce,
			"In Tungsten::update_heatcapacity():\nExtending heat capacity model outside temperature range! T < 300K");

		St.HeatCapacity = 24.943 - 7.72e4*pow(St.Temperature,-2) + 2.33e-3*St.Temperature + 1.18e-13*pow(St.Temperature,4);
	}else if( St.Temperature > 300 && St.Temperature < Ec.MeltingTemp ){ 
		// http://nvlpubs.nist.gov/nistpubs/jres/75a/jresv75an4p283_a1b.pdf
		St.HeatCapacity = 24.943 - 7.72e4*pow(St.Temperature,-2) + 2.33e-3*St.Temperature + 1.18e-13*pow(St.Temperature,4);
	}else{
		// http://webbook.nist.gov/cgi/inchi?ID=C7440337&Mask=2

		double t = St.Temperature/1000;
                St.HeatCapacity = 35.56404 -1.551741e-7*t + 2.915253e-8*pow(t,2) -1.891725e-9*pow(t,3)-4.107702e-7*pow(t,-2);
	}

//	E1_Debug("\n\nTemperature is : " << St.Temperature << "\nSt.Gas = " << St.Gas << "\nSt.Liquid = " << St.Liquid 
//				<< "\nCv of Solid: " << St.HeatCapacity/Ec.AtomicMass << "[kJ/(kg K)]"; );
	St.HeatCapacity = (St.HeatCapacity /(1000 * Ec.AtomicMass)); // Conversion kJ/(mol K) to kJ/( kg K ), AtomicMass [kg mol^-1]
}

void Tungsten::update_radius(){
	E_Debug("\tIn Tungsten::update_radius()\n\n");
	if( St.Temperature > 173 && St.Temperature <= 1500 ){
		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Tungsten::update_radius():\nExtending model outside temperature range!(from 738K to 1500K)");
		St.LinearExpansion = 1+(4.28*St.Temperature)*1e-6;
	}else if( St.Temperature > 1500 && St.Temperature < Ec.MeltingTemp ){
		St.LinearExpansion = 1+3.9003e-4+1.3896*1e-3-8.2797*1e-7*St.Temperature + 4.0557*1e-9*pow(St.Temperature,2)
					-1.2164*1e-12*pow(St.Temperature,3)+1.7034*1e-16*pow(St.Temperature,4);
		// 1.02648269488
	}else if( St.Temperature == Ec.MeltingTemp ){ // Model while it's melting
		// http://nvlpubs.nist.gov/nistpubs/jres/75a/jresv75an4p283_a1b.pdf
		St.LinearExpansion = 1.02105 + 0.03152*St.FusionEnergy/(Ec.LatentFusion*St.Mass);
		// 1.02648269488 - 1.05257238281 = 0.02608968793
	}else if( St.Temperature > Ec.MeltingTemp && St.Temperature <= St.SuperBoilingTemp){
		St.LinearExpansion = pow(1.18+6.20*1e-5*(St.Temperature-3680)+3.23*1e-8*pow((St.Temperature-3680),2),(1./3.));
	}
	St.Radius=St.UnheatedRadius*St.LinearExpansion;
	E1_Debug("\nTemperature = " << St.Temperature << "\n\nSt.LinearExpansion = " << St.LinearExpansion
			<< "\nSt.Radius = " << St.Radius);
	assert(St.Radius>0); // Assert radius is positive	
}

void Tungsten::update_vapourpressure(){
	E_Debug("\tIn Tungsten::update_vapourpressure()\n\n");
//	double AmbientPressure = 0;
	St.VapourPressure = 101325*pow(10,2.945 - 44094/St.Temperature + 1.3677*log10(St.Temperature)); // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
	
	// Answer is in pascals
}

