#ifndef __THREEVECTOR_H_INCLUDED__   // if plasmagrid.h hasn't been included yet...
#define __THREEVECTOR_H_INCLUDED__

#include <math.h>
#include <fstream>
#include <iostream>

class threevector{
	protected:
		double xcoord, ycoord, zcoord;

	public:
	
		threevector(); // Default constructor	
		threevector(double x, double y, double z); // Cartesian constructor
		threevector(double r, double theta, double phiorz, char type); // Polar constructor
	
		void print(); // Method to print out contents to screen	
		void print(std::ofstream &fout); // Method to print out contents to file
	
		// Modifier method to change a component of a threevector
		void setx(double xvalue); 
		void sety(double yvalue);
		void setz(double zvalue);
	
		// Access method for coordinates
		double getx()const;
		double gety()const;
		double getz()const;
	
		double square()const; 	// Calculating the square of a threevector	
		double mag3()const; 		// Calculating the magnitude 
			
		double gettheta()const; 	// Calculate theta
		double getphi()const; 		// Calculate phi in range 0 to 2pi
		
		threevector getunit()const; 	// Calculate a unit threevector
		
		threevector operator+(const threevector v_old)const; 	// Overload the + operator
		threevector operator-(const threevector v_old)const; 	// Overload the - operator
	
		void operator+=(threevector v_old);			// Overload the += operator
		void operator-=(threevector v_old);			// Overload the -= operator
		void operator=(threevector v_old); 			// Overload the = operator
		threevector operator*(double scalar)const; 		// Overload the * operator to multiply by a scalar
		double operator*(threevector v_old)const; 		// The dot product	
		threevector operator^(threevector v_old)const; 		// The threevector product
		
		friend std::ostream& operator<<(std::ostream& os, const threevector& vec);
};

threevector operator*(double scalar, const threevector& y);
threevector operator*(const bool scalar, const threevector& y);

#endif
