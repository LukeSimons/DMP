// 	Following work from:
//	Generating equally weighted test particles from the one-way flux of a drifting Maxwellian
//	T Makkonen, M I Airila & T Kurki-Suonio
//	http://iopscience.iop.org/article/10.1088/0031-8949/90/1/015204/data

// We need these headers
#include <cmath>
#include <cstdlib>
#include <random>	// for std::normal_distribution<> etc.

// forward definitions
double rand_mwts(double u, double sigma, std::mt19937 &mt);
double rand_mwts(double d, std::mt19937 &mt);
double solve_cubic(double d);
double rejection_function(double d, double t);

// Generates a number from the two parameter distribution
double rand_mwts(double u, double sigma, std::mt19937 &mt) {
  return rand_mwts(u/sigma,mt)*sigma;
}

// Generates a random number from the simplified distribution
double rand_mwts(double d, std::mt19937 &mt) {
  // Calculate the normalization for this d
  double maxval = rejection_function(d, solve_cubic(d));
  std::uniform_real_distribution<double> rad(0, 1); // Random number between 0 and 1
  double r1, r2;
  do {
    r1 = rad(mt);
    r2 = maxval*rad(mt);
  } while(r2 > rejection_function(d, r1));
  
  return r1 / (1.0 - r1);
}

// Solves the cubic equation used to find the maximum of the rejectin sampling function
double solve_cubic(double d) {

  // The variables below will be optimized away by the compiler
  // cubic equation
  double ca = 2.0;
  double cb = -4.0 - d;
  double cc = d;
  double cd = 1.0;
  // depressed cubic quation
  double p = (3.0*ca*cc - cb*cb)/(3.0*ca*ca);
  double q = (2.0*cb*cb*cb - 9.0*ca*cb*cc + 27*ca*ca*cd) / (27*ca*ca*ca);


  // calculate the value
  double sqrt_term = sqrt(-p/3.0);
  double root = 2.0*sqrt_term*cos(acos(3.0*q/(2.0*p*sqrt_term))/3.0 - 2.0*M_PI/3.0); 
  return root - cb/(3.0*ca);
}

// Gives the unnormalized function used in the rejection sampling algorithm
double rejection_function(double d, double t) {
  double s = 1.0-t;
  return t*exp(-(t/s - d)*(t/s - d)/2.0)/(s*s*s);
}

