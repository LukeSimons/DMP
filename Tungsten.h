#ifndef __TUNGSTEN_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __TUNGSTEN_H_INCLUDED__

//#include <iostream>
#include "Matter.h"

class Tungsten: public Matter{

	private:

		// Functions called by Tungsten::update()
		void update_radius		();
		void update_heatcapacity 	();
		void update_vapourpressure	();

	public:
		// Constructors
		Tungsten();
		Tungsten(double radius);
		Tungsten(double radius, double tempin);
		Tungsten(double radius, double tempin, std::array<char,4> &constmodels);

		// Destructor
		~Tungsten(){};
		
		// Change Properties; Mass and Temperature
		void tungsten_defaults		();
};

#endif
