
// RANDom number from "MaxWellian Trough Surface" distribution
// Toni Makkonen 17.1.2013
// toni.makkonen@aalto.fi

// Generates random number from the two parameter distribution
// u is the flow velocity
// sigma is the thermal velocity (square root of variance of velocity distribution)
// mt is the Mersenne Twister random number generator
double rand_mwts(double u, double sigma, std::mt19937 &mt);

// Generates random number from the simplified one parameter distribution in the paper
// mt is the Mersenne Twister random number generator
double rand_mwts(double d, std::mt19937 &mt);
