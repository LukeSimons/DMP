#ifndef __MATTER_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __MATTER_H_INCLUDED__

#include "GrainStructs.h"// Contains the structure element for storing properties of material
#include "Constants.h"	// Contains general physical constants
#include "Functions.h"  // sec(Te,'f') function used by HeatingModel.cpp
#include "threevector.h"

#include <iostream>	// I/O operations and debugging
#include <fstream> 	// Open emissivity filestream
#include <sstream>	// Convert temp to file number
#include <assert.h>	// Assertion errors
#include <math.h> 	// Round
#include <array>	// std::array
#include <limits>	// std::numeric_limits<double>::max()
#include <string.h>	// strchr("",std::string)

class Matter{

       	// Functions used by HeatingMatter.cpp
	protected:
		// Constructors
		Matter(const ElementConsts *elementconsts);
		Matter(double rad, const ElementConsts *elementconsts);
		Matter(double rad, double temp, const ElementConsts *elementconsts);
		Matter(double rad, double temp, const ElementConsts *elementconsts, std::array <char,4> &constmodels);
 
		struct GrainData 		St;
		const struct ElementConsts	Ec;
		double PreBoilMass;			// A variable needed to remember the mass before boiling...
		std::array<char,4> ConstModels;		// Constant Models variation with Temperature turned on of possibly 4

		// Functions called by Matter::update(). These are element dependant and must be defined in the child class.
		virtual void update_radius		()=0;
		virtual void update_heatcapacity 	()=0;
		virtual void update_vapourpressure	()=0;

		// Update geometric properties of matter:
		// Volume, Surface area and Density. Also initialises mass
		void update_dim();
		void update_emissivity();
		void update_state(double EnergyIn);
		void update_boilingtemp();

	public:
		// Destructor
		virtual ~Matter			(){M_Debug("M_Debug was on.\n");};

		// Change geometric and variable properties of matter; 
		// volume, surface area, density, heat capacity, emissivity and vapour pressure
		// Think about putting 'update()' in HeatingMatter.cpp
		void update();
		void update_models(char emissivmodel, char linexpanmodel, char heatcapacity, char boilingmodel);
		void update_models(std::array<char,4> &constmodels);
		void update_mass(double LostMass); // Mass lost in Kilogrammes
		// Takes argument of amount of energy lost in Kilo Joules and changes temperature in Degrees
		void update_temperature(double EnergyIn);

		void update_motion(threevector &changeinposition, threevector &changeinvelocity, double Rotation);
		// Resolve region of rapid charge variation
		void update_charge(double charge, double potential, double deltas, double deltat);

		// Getter methods
		char get_elem			()const{ return Ec.Elem;			};
		double get_meltingtemp		()const{ return Ec.MeltingTemp;			};
		double get_boilingtemp		()const{ return Ec.BoilingTemp;			};
		double get_bondenergy		()const{ return Ec.BondEnergy;			};
		double get_workfunction		()const{ return Ec.WorkFunction;		};
		double get_heattransair		()const{ return Ec.HeatTransAir;		};
		double get_atomicmass		()const{ return Ec.AtomicMass;			};
		double get_surfacetension	()const{ return Ec.SurfaceTension;		};

		bool is_gas			()const{ return St.Gas;				};
		bool is_liquid			()const{ return St.Liquid;			};
		bool is_split			()const{ return St.Breakup;			};
		bool is_positive		()const{ return St.Positive;			};
		bool get_c		   (int i)const{ assert(i < 4); return ConstModels[i];	};
		double get_superboilingtemp	()const{ return St.SuperBoilingTemp;		};
		double get_mass			()const{ return St.Mass;			};
		double get_density		()const{ return St.Density;			};
		double get_radius		()const{ return St.Radius;			};
		double get_heatcapacity		()const{ return St.HeatCapacity;		};
		double get_temperature		()const{ return St.Temperature; 		};
		double get_fusionenergy		()const{ return St.FusionEnergy;		};
	        double get_vapourenergy         ()const{ return St.VapourEnergy;		};
		double get_latentvapour		()const{ return Ec.LatentVapour;		};
		double get_latentfusion		()const{ return Ec.LatentFusion;		};
		double get_emissivity		()const{ return St.Emissivity; 			};
		double get_surfacearea		()const{ return St.SurfaceArea;			};
		double get_vapourpressure	()const{ return St.VapourPressure;		};
		double get_linearexpansion	()const{ return St.LinearExpansion;		};
		double get_deltasec		()const{ return St.DeltaSec;			};
		double get_deltatherm		()const{ return St.DeltaTherm;			};
		double get_deltatot		()const{ return (St.DeltaTherm + St.DeltaSec);	};
		double get_potential		()const{ return St.Potential;			};
		double get_rotationalfreq	()const{ return St.RotationalFrequency;		};
		threevector get_velocity	()const{ return St.DustVelocity;		};
		threevector get_position	()const{ return St.DustPosition;		};


		void set_potential		(double potential){ St.Potential = potential;	};
};

#endif
