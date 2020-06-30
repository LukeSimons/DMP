/** @file rand_mwts.h
 *  @brief  RANDom number from "MaxWellian Trough Surface" distribution
 *  See "Generating equally weighted test particles from the one-way flux of 
 *  a drifting Maxwellian". T Makkonen, M I Airila & T Kurki-Suonio
 *  http://iopscience.iop.org/article/10.1088/0031-8949/90/1/015204/data
 *  
 *  @author Toni Makkonen (toni.makkonen@aalto.fi)
 *  @date 17.1.2013
 *  @bug No known bugs
 */

//!< Generates random number from the two parameter distribution
//!< u is the flow velocity
//!< sigma is the thermal velocity (square root of variance of velocity 
//!< distribution)
//!< mt is the Mersenne Twister random number generator
double rand_mwts(double u, double sigma, std::mt19937 &mt);

//!< Generates random number from the simplified one parameter distribution in 
//!< the paper mt is the Mersenne Twister random number generator
double rand_mwts(double d, std::mt19937 &mt);
