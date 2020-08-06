/** @file main.cpp
 *  @brief DiMPl program to calculate behaviour of conducting spheres in
 *  plasmas with background flow and magnetic fields
 *
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs
 */

#include <omp.h>     //!< For omp parallelisation
#include <iostream>  //!< for std::cout
#include <array>
#include <random>    //!< for std::normal_distribution<> etc.
#include <fstream>   //!< for std::ofstream
#include <ctime>     //!< for clock()
#include <math.h>    //!< for fabs()
#include <sstream>   //!< for std::stringstream
#include <assert.h>  //!< for assert()
#include <algorithm> //!< for std::min()
#include <vector> //!< To store text file of unknown length

#include "Switches.h"   //!< File containing switches defining DiMPl behaviour
#include "Constants.h"  //!< Define Pre-processor Directives and constants
#include "threevector.h"//!< For threevector class
#include "rand_mwts.h"  //!< Functions to generate correct 1-way maxwellian flux
#include "Field_Point.h" //!< For Field_Point class (contains value and position vector)
#include "Field_Map.h" //!< For Field_Map class (contains details of map of field point objects)

using namespace dimplconsts;

/** @brief display a help message to the user on command line options
 *  @param name is the text passed by user to the function
 *
 *  Display to the user what options are available for running the program
 */
static void show_usage(std::string name){
    std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
    << "\n\nOptions:\n"
    << "\t-h,--help\t\t\tShow this help message\n\n"
    << "\t-r,--radius RADIUS\t\t(m), Specify radius of Dust grain\n"
    << "\t\tRadius(=1e-6m) DEFAULT,\t\tBy Default, simulate sphere of size 1um "
        << "in radius\n\n"
    << "\t-s,--spin SPIN\t\t(hz), Specify the initial rotation frequency of "
        << "Dust grain alligned with magnetic field axis\n"
    << "\t\tSpin(=0.0hz) DEFAULT,\t\tBy Default, simulate sphere with "
        << "initially zero z angular momentum\n\n"
    << "\t-a1,--semix SEMIX\t\t(arb), Specify the semi-axis for x in dust "
        << "radii\n"
    << "\t\ta1(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
    << "\t-a2,--semiy SEMIY\t\t(arb), Specify the semi-axis for y in dust "
        << "radii\n"
    << "\t\ta2(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
    << "\t-a3,--semiz SEMIZ\t\t(arb), Specify the semi-axis for z in dust "
        << "radii\n"
    << "\t\ta3(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
    << "\t-d,--density DENSITY\t\t(kgm^-^3), Specify density of Dust grain\n"
    << "\t\tDensity(=19600kgm^-^3) DEFAULT,\tBy Default, Tungsten Density\n\n"
    << "\t-ex,--efieldx EFIELDX\t(double), Specify the x component of a "
        << "background electric field\n"
    << "\t-ey,--efieldx EFIELDY\t(double), Specify the y component of a "
        << "background electric field\n"
    << "\t-ez,--efieldx EFIELDZ\t(double), Specify the z component of a "
        << "background electric field\n\n"
    << "\t-p,--potential POTENTIAL\t(double), Specify the potential of Dust "
        << "grain normalised to electron temperature\n"
    << "\t\tPotential(=-2.5eV) DEFAULT,\tBy Default, OML Potential in Ti=Te "
        << "Hydrogen plasma\n\n"
    << "\t-m,--magfield MAGFIELD\t\t(T), Specify the magnetic field "
        << "(z direction)\n"
    << "\t\tBMag(=1.0T) DEFAULT,\t\tBy Default, magnetic field is 1.0T upwards "
        << "in vertical z\n\n"
    << "\t-n,--normalised NORMALISED\t(bool), whether normalisation "
        << "(following Sonmor & Laframboise) is on or off\n"
    << "\t\tNormVars(=0) DEFAULT,\tFalse, By Default use Tesla and "
        << "electron volts\n\n"
    << "\t-te,--etemp ETEMP\t\t(eV), Specify the temperature of plasma "
        << "electrons\n"
    << "\t\teTemp(=1eV) DEFAULT,\n\n"
    << "\t-ne,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron "
        << "density\n"
    << "\t\teDensity(=1e18m^-^3) DEFAULT,\tTokamak density\n\n"
    << "\t-mi,--imass IMASS\t\t(u), Specify the mass of Ions in Daltons\n"
    << "\t\tiMass(=1u) DEFAULT,\n\n"
    << "\t-ti,--itemp ITEMP\t\t(eV), Specify the temperature of Ions\n"
    << "\t\tiTemp(=1eV) DEFAULT,\n\n"
    << "\t-ni,--idensity IDENSITY\t\t(m^-^3), Specify the plasma ion density\n"
    << "\t\tiDensity(=1e18m^-^3) DEFAULT,\tTokamak density\n\n"
    << "\t-c,--ichance ICHANCE\t\t\t(eV), Specify the probability of "
        << "generating an Ion\n"
    << "\t\tiChance(=-0.5) DEFAULT,\tFicticious ion generation probability: "
        << "i.e Self-consistently generate ions & electrons\n\n"
    << "\t-u,--zmaxcoeff ZMAXCOEFF\t(double), The upper limit of simulation "
        << "domain as number of Coulomb Interaction lengths\n"
    << "\t\tzMaxCoeff(=1.0) DEFAULT,\tDistance from "
        << "(0,0,Radius) plane to max of simulation domain in dust radii\n\n"
    << "\t-l,--zmincoeff ZMINCOEFF\t(double), The lower limit of simulation "
        << "domain as number of Coulomb Interaction lengths\n"
    << "\t\tzMaxCoeff(=1.0) DEFAULT,\tDistance from "
        << "(0,0,Radius) plane to min of simulation domain in dust radii\n\n"
    << "\t-z,--zboundforce ZBOUNDFORCE\t(double), Force the absolute value of "
        << "simulation domain upper and lower boundaries\n"
    << "\t\tZBoundForce(=0.0) DEFAULT,\tBy Default, use -u and -l to determine "
        << "height of injection plane\n\n"
    << "\t-b,--impactpar IMPACTPAR\t(double), Specify the radial limit of "
        << "simulation domain as multiple of gyro-radii and debye lengths\n"
    << "\t\tImpactPar(=2.0) DEFAULT,\tBy Default, Radial extent of injection "
        << "is two gyro-radii plus two debye lengths from dust\n\n"
    << "\t-f,--forceimppar FORCEIMPPAR\t(double), Force the absolute value of "
        << "simulation radial distance\n"
    << "\t\tForceImpPar(=0.0) DEFAULT,\tBy Default, use -b to determine radial "
        << "extent of injection plane\n\n"
    << "\t-i,--imax IMAX\t\t\t(int), Specify the number of particles to be "
        << "launched\n"
    << "\t\timax(=10000) DEFAULT,\t\tBy Default, inject 10,000 particles\n\n"
    << "\t-j,--jmax JMAX\t\t\t(int), Specify the number of particles to be "
        << "collected (not exceeding imax)\n"
    << "\t\tjmax(=5000) DEFAULT,\t\tBy Default, stop simulation if 5,000 "
        << "particles are collected\n\n"
    << "\t-rm,--rmax RMAX\t\t\t(int), Specify the number of reflections in z "
        << "axis before assuming orbit misses\n"
    << "\t\trmax(=15) DEFAULT,\t\tBy Default, stop simulation if Particle "
        << "reflects in z axis 15 times\n\n"
    << "\t-t,--timestep TIMESTEP\t\t(double), Specify the multiplicative "
        << "factor for the time step\n"
    << "\t\tTimeStepFactor(=0.0005) DEFAULT,\n\n"
    << "\t-no,--number NUMBER\t\t\t(int), Specify the number of particles to "
        << "be captured before saving\n"
    << "\t\tNumber(=1000) DEFAULT,\t\tBy Default, store average values after "
        << "collecting 1000 particles\n\n"
    << "\t-v,--driftvel DRIFTVEL\t\t(m s^-^1), Specify the drift velocity of "
        << "the plasma\n"
    << "\t\tDriftVel(=0.0) DEFAULT,\t\tBy Default, No drift velocity\n\n"
    << "\t-se,--seed SEED\t\t\t(double), Specify the seed for the random "
        << "number generator\n"
    << "\t\tSeed(=1.0) DEFAULT,\t\tBy Default, Seed is 1.0\n\n"
    << "\t-sa,--saves SAVES\t\t\t(int), Specify the number of particles before "
        << "saving a run\n"
    << "\t\tSaves(=1000) DEFAULT,\t\tBy Default, Save data to Meta datafile "
        << "once per 1000.0 particles\n\n"
    << "\t-o,--output OUTPUT\t\t(string), Specify the suffix of the output "
        << "file\n"
    << "\t\tsuffix(='.txt') DEFAULT,\tBy Default, Save data to "
        << "Data/DiMPl.txt\n\n"
    << std::endl;
    #if defined CUSTOM_POTENTIAL
    std::cerr <<"\n\nAdditional Options from CUSTOM_POTENTIAL.\n"
    << "\t-cf, --customfield CUSTOMFIELD\t\t(string), Specify the Custom "
        << "Field file name to be used\n"
    << "\t\tcustom electric potential file name"
        << "(='Custom_Fields/Default_Field.txt') DEFAULT, \tBy Default, Load "
	<< "file Custom_Fields/Default_Field.txt\n\n"
    << std::endl;
    #endif
    #if defined VARIABLE_CSCALE
    std::cerr << "\n\nAdditional Options from VARIABLE_CSCALE!\n"
    << "\t-cs,--chargescale CHARGESCALE\t\t\t((double), Specify the scale of "
        << "the charging\n"
    << "\t\tChargeScale(=0) DEFAULT,\t\tBy Default, charging scale is "
        << "calculated from Potential~1\n\n"
    << std::endl;
    #endif
    #if defined VARIABLE_CSCALE || defined VARIABLE_ASCALE
    std::cerr << "\n\nAdditional Options from VARIABLE_CSCALE OR "
        << "VARIABLE_ASCALE!\n"
    << "\t-jm,--jmin JMIN\t\t\t(int), Specify the number of particles to be "
        << "collected before dynamic saving\n"
    << "\t\tjmin(=20) DEFAULT,\t\tBy Default, re-assess dynamic scale after "
        << "collecting 20 particles\n\n"
    << "\t-jf,--jfin JMIN\t\t\t(int), Specify the number of particles to be "
        << "collected in final save\n"
    << "\t\tjfin(=100) DEFAULT,\t\tBy Default, collect 100 particles in the "
        << "final save\n\n"
    << std::endl;
    #endif

}

/** @brief Function to process user command line input
 *  @param argc The number of arguments input by the user
 *  @param argv The character array which contains the users input
 *  @param i A counter indicating the argument which is being processed
 *  @param ss0 The string stream which is used to process command line input
 *  @param Temp The variable of type T which is being filled with user input
 *  @return An integer representing status of input retreival
 *
 *  Process user command line input and return status of success of failure
 */
template<typename T> int InputFunction(int &argc, char* argv[], int &i,
std::stringstream &ss0, T &Temp){
    if (i + 1 < argc) { //!< Make sure we aren't at the end of argv!
        i+=1; //!< Increment 'i' so we don't get the next argv[i].
        ss0 << argv[i];
        ss0 >> Temp;
        ss0.clear(); ss0.str("");
        return 0;
    }else{ //!< Uh-oh, there was no argument to the destination option.
        std::cerr << "\noption requires argument." << std::endl;
        return 1;
    }

}

#ifdef COLLISIONS
/** @brief function to calculate the probability of a collision
 *  @param time the time period over which to calculate probability
 *  @param velocity of incident particle
 *  @param meanfreepath scale length defining collision probability
 *  @param mt mersenne twister which generates random numbers
 *  @return A boolean true if particle collides and false otherwise
 *
 *  Function to determine whether a particle undergoes a charge exchange
 *  collision with a neutral particle.
 */
bool collision_probability(double time, double velocity, double meanfreepath,
    std::mt19937 &mt)
{
    //!< Random uniform Distribution
    std::uniform_real_distribution<double> rad(0.0, 1.0);
    double u = rad(mt); //!< Random number between 0 and 1
    bool returnvalue = (u < 1.0-exp(-time*velocity/meanfreepath));
    return returnvalue;
}
#endif

#ifdef SPHERICAL_INJECTION
/** @brief convert a velocity from spherical to cartesian coordinates
 *  @param v pointer to the three vector of velocity in spherical coordinates
 *  @param phi the angle phi in the spherical coordinate system
 *  @param theta the angle theta in the spherical coordinate system
 *  @param vx pointer to x component of cartesian velocity
 *  @param vy pointer to y component of cartesian velocity
 *  @param vz pointer to z component of cartesian velocity
 *  @author d. m. thomas (drew.thomas07@imperial.ac.uk)
 *
 *  convert a velocity in spherical coordinates to a velocity in cartesian
 *  coordinates, at a position given in spherical coordinates.
 */
void velocity_in_cartesian_coords(double* v, double phi, double theta,
    double* vx, double* vy, double* vz)
{
    *vx = (v[0] * cos(phi) * sin(theta))
          - (v[1] * sin(phi)) + (v[2] * cos(phi) * cos(theta));
    *vy = (v[0] * sin(phi) * sin(theta))
          + (v[1] * cos(phi)) + (v[2] * sin(phi) * cos(theta));
    *vz = (v[0] * cos(theta)) - (v[2] * sin(theta));
}

/** @brief Generate pair of correlated normal random variates
 *  @param sd standard deviation of the distribution
 *  @param nrv1 pointer to the first random variable
 *  @param nrv2 pointer to the second random variable
 *  @param mt mersenne twister which generates random numbers
 *  @author D. M. Thomas (drew.thomas07@imperial.ac.uk)
 *
 *  Generate a pair of normal random variates with mean zero and store them
 *  in `nrv1` and `nrv2`.
 */
void nrv_pair(double sd, double* nrv1, double* nrv2, std::mt19937 &mt)
{
    double r;
    double theta;
    long double urv1;

    std::uniform_real_distribution<double> rad(0.0, 1.0); // IONS

    long double urv2 = rad(mt);
    //!< Calculate the variates with the Box-Muller transform's polar form.
    do {
        urv1 = rad(mt);
    } while (urv1 == 0.0);  //!< `rad` very occasionally returns zero!
    r = sqrt(-2.0 * (double) logl(urv1));
    theta = 2.0 * PI * urv2;
    *nrv1 = sd * r * cos(theta);
    *nrv2 = sd * r * sin(theta);
}

/** @brief probability distribution function of inclination angle theta, under a
 *  flow mach u.
 *  @param theta the angle of inclination of the distribution
 *  @param u flow velocity normalised to thermal velocity
 *  @param sigma standard deviation of the distribution
 *  @return double
 *
 *  probability distribution from inclination angle \p theta, under a flow mach
 *  \p u
 */
double thetaPDF(double theta, double u, double sigma)
{
    return sin(theta)*(1.0+erf(-1*u/(tan(theta)*sigma
        *sqrt(2.0+2/(tan(theta)*tan(theta))))));
}

/** @brief find the maximum value of the theta probability distribution
 *  function
 *  @param u flow velocity normalised to thermal velocity
 *  @param sigma standard deviation of the distribution
 *  @param tol the tolerance of maximisation function, default=1e-6
 *  probability distribution from inclination angle \p theta, under a flow mach
 *  \p u
 *
 *  Evaluate maximum value of theta pdf, for use in rejection sampling
 */
double thetaPDFMax(double u, double sigma, double tol=1e-6)
{
    double d_theta = 0.1;
    double theta_0 = 0.0;
    double theta_1 = PI;
    double theta = 0.0;
    while (d_theta > tol)
    {
        theta = theta_0;
        while (theta < theta_1)
        {
            if (thetaPDF(theta,u,sigma) >= thetaPDF(theta-d_theta,u,sigma)
                && (thetaPDF(theta,u,sigma) >= thetaPDF(theta+d_theta,u,sigma)))
            {
                theta_0 = theta-d_theta;
                theta_1 = theta+d_theta;
                d_theta /= 10.0;
                break;
            }else{
                theta += d_theta;
            }
        }
    }
    return thetaPDF(theta_1-10.0*d_theta,u,sigma);
}


#endif

/** @brief Generate randomly distributed initial positions and velocities
 *  @param Position reference to three vector with generated position
 *  @param Velocity reference to three vector with generated velocity
 *  @param ImpactParameter radius of injection plane (circle or sphere)
 *  @param ProbUpper probability of injecting particles on upper plane
 *  @param zmin distance to lower injection plane
 *  @param zmax distance to upper injection plane
 *  @param DriftNorm the background flow velocity, shifted distribution
 *  @param ThermalVel the thermal velocity defining width of distribution
 *  @param mt mersenne twister which generates random numbers
 *
 *  Generate the random initial positions and velocities of particles over
 *  distributions determined by switches. Default is uniformly distributed
 *  positions over circular area with velocities given by rand_mwts().
 *  Spherical injection and point injection are also possible
 */
void GenerateOrbit(threevector &Position, threevector &Velocity,
    const double &ImpactParameter, const double &ProbUpper, const double &zmin,
    const double zmax, const double DriftNorm, const double ThermalVel,
    const double xThetaPDFMax, std::mt19937 &mt){

    // ***** DEFINE RANDOM NUMBER GENERATOR ***** //
    //!< See
    //!< https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    //!< search for
    //!< "Each component of the velocity vector has a normal distribution"

    std::uniform_real_distribution<double> rad(0.0, 1.0);
    #ifdef SPHERICAL_INJECTION
        // ***** RANDOMISE POSITION SPHERICALLY ***** //
        double phi_pos = 2.0*PI*rad(mt);
        double theta_pos;
        bool accepted = false;

        // rejection sampling to select theta coordinate
        while (!accepted)
        {
            theta_pos = PI*rad(mt);
            double x = xThetaPDFMax*rad(mt);
            if (x < thetaPDF(theta_pos,DriftNorm,ThermalVel)){
                accepted = true;
            }
        }

        Position.setx(ImpactParameter*sin(theta_pos)*cos(phi_pos));
        Position.sety(ImpactParameter*sin(theta_pos)*sin(phi_pos));
        Position.setz(ImpactParameter*cos(theta_pos));


        // ***** RANDOMISE VELOCITY SPHERICALLY ***** //
/*  POT Style velocity generation
        double v[3],v_temp[3];  //!< 0 is v_r, 1 is v_phi, and 2 is v_theta
        nrv_pair(ThermalVel, &(v[0]), &(v[1]),mt);
        nrv_pair(ThermalVel, &(v[2]), &(v[2]),mt);

        if (v[0] > 0.0) {
            //!< No point injecting a particle at the outer boundary
            //!<   with positive radial velocity.
            v[0] = -v[0];
        }
*/

        double invel = rand_mwts(0.0,ThermalVel,mt);
        std::normal_distribution<double> Gaussdist(0.0,ThermalVel);
        double v[3], v_temp[3];
        v[0]= -invel;
        v[1]= Gaussdist(mt);
        v[2]= Gaussdist(mt);

        /* Translate the injection velocity in spherical coordinates into
           the desired Cartesian coordinates. Luckily I already had this
           worked out for John and me's spherical absorber problem! */
        velocity_in_cartesian_coords(v, phi_pos, theta_pos,
            &(v_temp[0]), &(v_temp[1]), &(v_temp[2]));

        Velocity.setx(v_temp[0]);
        Velocity.sety(v_temp[1]);
        Velocity.setz(v_temp[2]+DriftNorm);

        //!< For case of flow, if the injected particle is
        //!< expected to leave the simulation domain immediately,
        //!< reflect the position
        if( Velocity.getz()*Position.getz() > 0.0 )
            Position.setz(-1.0*Position.getz());

    #elif defined(POINT_INJECTION)
        // ***** INJECT AT SINGLE POINT WITH DRIFT VELOCITY ***** //
        Position.setx(ImpactParameter);
        Position.sety(0.0);
        Position.setz(zmax);
        double phi_pos = 2.0*PI*rad(mt);
        double theta_pos = PI*rad(mt);
        if( DriftNorm == 0.0 ){
            Velocity.setx(ThermalVel*sin(theta_pos)*cos(phi_pos));
            Velocity.sety(ThermalVel*sin(theta_pos)*sin(phi_pos));
            Velocity.setz(-1.0*ThermalVel*fabs(cos(theta_pos)));
        }else{
            Velocity.setx(0.0);
            Velocity.sety(0.0);
            Velocity.setz(DriftNorm);
        }
    #else
        // ***** RANDOMISE POSITION CYLINDRICALLY ***** //
        double radial_pos = ImpactParameter*sqrt(rad(mt));
        double theta_pos = 2*PI*rad(mt);
        Position.setx(radial_pos*cos(theta_pos));
        Position.sety(radial_pos*sin(theta_pos));

        // ***** RANDOMISE VELOCITY CYLINDRICALLY ***** //
        // Following work from:
        // Generating equally weighted test particles from the one-way flux of
        // a drifting Maxwellian
        // T Makkonen, M I Airila & T Kurki-Suonio
        // http://iopscience.iop.org/article/10.1088/0031-8949/90/1/015204/data

        double invel(0.0);

        if( rad(mt) > ProbUpper ){
            Position.setz(zmax);
            // Generate z-velocity here to determine starting position
            invel = -rand_mwts(-DriftNorm,ThermalVel,mt);
        }else{
            Position.setz(zmin);
            // Generate z-velocity here to determine starting position
            invel = rand_mwts(DriftNorm,ThermalVel,mt);
        }
        std::normal_distribution<double> Gaussdist(0.0,ThermalVel);

        Velocity.setx(Gaussdist(mt));
        Velocity.sety(Gaussdist(mt));
        Velocity.setz(invel);
    #endif
}

/** @brief Implement Boris algorithm to solve particle trajectories
 *  @param MASS the mass of the particle
 *  @param Efield reference to the three vector of electric field
 *  @param BField reference to the three vector of magnetic field
 *  @param dt the time step used in solving algorithm
 *  @param Velocity reference to the three vector of particle velocity
 *  @param SPEC_CHARGE the charge of plasma particle species
 *
 *  Solve particle trajectories explicitly by the Boris algorithm. updates
 *  velocity using the Boris method, Birdsall, Plasma Physics via Computer
 *  Simulation, p.62
 */
static void UpdateVelocityBoris(double MASS, const threevector& Efield,
    const threevector& BField, double dt, threevector &Velocity,
    const int SPEC_CHARGE){

    /*magnitude of t, squared*/
    double t_mag2 = (BField*0.5*SPEC_CHARGE*dt).square();

    /*v prime*/
    threevector v_prime = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt
        + ((Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt)^BField
        *0.5*SPEC_CHARGE*dt);

    /*v prime*/
    threevector v_plus = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt
        + (v_prime^(2.0*BField*0.5*SPEC_CHARGE*dt*(1.0/(1.0+t_mag2))));

    /*v n+1/2*/
    Velocity = v_plus + Efield*0.5*(SPEC_CHARGE/MASS)*dt;
}

/** @brief Calculate the electric field due to central coulomb potential
 *  @param Position reference to the three vector of particle position
 *  @param Charge the charge of the particle
 *  @param COULOMB_NORM the normalisation of the electric field
 *  @return threevector of the coulomb electric field at \p Position
 *
 *  Calculate the coulomb electric field for a \p Charge at position
 *  \p Position with normalisation \p COULOMB_NORM
 */
#ifdef COULOMB_POTENTIAL
threevector CoulombField(const threevector &Position, double Charge,
    double COULOMB_NORM){
    return (COULOMB_NORM*Charge/Position.square())*Position.getunit();
}
#endif

/** @brief Calculate the electric field due to central debye-huckel potential
 *  @param Position reference to the three vector of particle position
 *  @param Charge the charge of the particle
 *  @param DebyeLength the debye length of the plasma
 *  @param COULOMB_NORM the normalisation of the electric field
 *  @return threevector of the debye-huckel electric field at \p Position
 *
 *  Calculate the debye-huckel electric field for a \p Charge at position
 *  \p Position with normalisation \p COULOMB_NORM and debye length
 *  \p DebyeLength
 */
#ifdef DEBYE_POTENTIAL
threevector DebyeHuckelField(const threevector &Position, double Charge,
    double DebyeLength, double COULOMB_NORM){
    if(Charge==0.0) return threevector(0.0,0.0,0.0);
    return COULOMB_NORM*Charge*(1.0/(Position.mag3()))
        *exp(-(Position.mag3()-1.0)/DebyeLength)
        *(1.0/Position.mag3()+1.0/DebyeLength)*Position.getunit();
}
#endif

/** @brief Calculate the electric field due to central parabolic potential
 *  @param Position reference to the three vector of particle position
 *  @param Charge the charge of the particle
 *  @param SheathLength the length scale of the parabola
 *  @param COULOMB_NORM the normalisation of the electric field
 *  @return threevector of the parabolic electric field at \p Position
 *
 *  Calculate the parabolic electric field for a \p Charge at position
 *  \p Position with normalisation \p COULOMB_NORM and length scale
 *  \p SheathLength
 */
#ifdef PARABOLIC_POTENTIAL
threevector ParabolicField(const threevector &Position, double Charge,
    double SheathLength, double COULOMB_NORM){
    if( (Position.mag3()-1.0) > SheathLength ) return threevector(0.0,0.0,0.0);
    return (-2.0*COULOMB_NORM*Charge*(Position.mag3()-1.0-SheathLength)
        /(SheathLength*SheathLength))*Position.getunit();
}
#endif

/** @brief Calculate the electric field from a custom potential field file
 *  @param Position reference to the three vector of particle position
 *  @return threevector of the custom electric field at \p Position
 *
 *  Find the approximate custom electric field at position
 *  \p Position
 */
#ifdef CUSTOM_POTENTIAL
threevector CustomField(const threevector &Position, Field_Map &this_field_map){
    return this_field_map.find_approx_value(Position, false);
}
#endif

/** @brief Main function defining DiMPl program
 *  @param argc the number of command line inputs given by the user
 *  @param argv the character array defining user command line inputs
 *  @return return status defining the status of the program
 *
 *  DiMPl calculates the trajectories of ions and electrons in the vicinity
 *  of a charged conducting sphere. See README and docs for more information.
 */
int main(int argc, char* argv[]){

    // ***** TIMER AND FILE DECLERATIONS        ***** //
    clock_t begin = clock();
    std::string filename = "Data/DiMPl";
    std::string suffix  = ".txt";
    DECLARE_TRACK();  //!< Data file for recording all positional information
    DECLARE_CHARGE(); //!< Data file for recording dynamic charging behaviour
    DECLARE_MOM();    //!< Data file for momentum of missed particles
    DECLARE_AVEL();   //!< Data file for collected angular momentum
    DECLARE_LMOM();   //!< Data file for collected linear momentum
    DECLARE_CHA();    //!< Data file for charge
    DECLARE_SPOS();   //!< Data file for starting positions
    DECLARE_EPOS();   //!< Data file for end positions
    DECLARE_APP();    //!< Date file for closest approaches
    DECLARE_CURR();   //!< Data file for the Currents
    DECLARE_TOT();    //!< Data file for the particle totals
    DECLARE_REF();    //!< Data file for the reflections
    DECLARE_AXIS();   //!< Data file for tracking free dust axis
    std::ofstream RunDataFile;   //!< Data file for containing the run
    std::ofstream InputDataFile; //!< Data file for containing the input
    #ifdef CUSTOM_POTENTIAL
    const std::string FIELD_DIRECTORY = "Custom_Fields";
    const std::string DIRECTORY_DELIM = "/";
    const std::string DEFAULT_FIELD_FILENAME = "Default_Field.txt";
    std::string input_custom_filename = DEFAULT_FIELD_FILENAME;
    #endif

    // ************************************************** //


    // ***** DEFINE DUST PARAMETERS         ***** //
    double Radius    = 1e-6;  //!< m, Radius of dust
    double Spin      = 0.0;   //!< hz, Initial rotation rate
    double Density   = 19600; //!< kg m^-^3, Tungsten
    double Potential = -2.5;  //!< Normalised Potential,
    double BMag      = 1.0;   //!< Tesla, Magnitude of magnetic field
    double BMagIn    = BMag;  //!< (arb), input magnetic field, in Tesla
    bool   NormVars  = false; //!< Normalised to Sonmor & Laframboise?
    double a1        = 1.0;   //!< Semi-axis for x in dust-grain radii
    double a2        = 1.0;   //!< Semi-axis for y in dust-grain radii
    double a3        = 1.0;   //!< Semi-axis for z in dust-grain radii
    #ifdef FREE_AXIS
        threevector DustAxis(0.0,0.0,0.0);
    #endif

    // ************************************************** //


    // ***** DEFINE PLASMA PARAMETERS       ***** //
    double EFieldx      = 0.0;  //!< V/m, x component of electric field
    double EFieldy      = 0.0;  //!< V/m, x component of electric field
    double EFieldz      = 0.0;  //!< V/m, x component of electric field
    double iMass        = 1.00784;  //!< Ion Mass, u
    double iTemp        = 1.0;  //!< Ion Temperature, eV
    double eTemp        = 1.0;  //!< Electron Temperature, eV
    double eDensity     = 1e18; //!< m^(-3), Electron density
    double iDensity     = 1e18; //!< m^(-3), Ion density
    double DriftVel     = 0.0;  //!< m s^-1, This is the Temperature of the species being considered!<
    double zMaxCoeff    = 1.0;  //!< Arb, Number of Interaction distances from 0,0 plane to max of simulation domain
    double zMinCoeff    = 1.0;  //!< Arb, Number of Interaction distances from 0,0 plane to min of simulation domain
    double ZBoundForce  = 0.0;  //!< Arb, Number of dust grain radii to vertical edge of simulation domain
    double ImpactPar    = 2.0;  //!< Species gyro-radii, Multiplicative factor for the Impact Parameter
    double ForceImpPar  = 0.0;  //!< Number of dust grain radii to radial edge of simulation domain
    double iChance      = -0.5; //!< Manually set probability of Generating an ion, negative will cause self-consistent
    unsigned long long imax = 10000;//!< Arb, Maximum number of particles to be launched
    unsigned long long jmax = 5000; //!< Arb, Number of particles to be collected
    #if defined VARIABLE_CSCALE || defined VARIABLE_ASCALE
    unsigned long long jmin = 20;    //!< Arb, Min particles before averaging
    unsigned long long jfin = 100;   //!< Arb, Min particles in last save
    #endif
    #if defined VARIABLE_CSCALE
    double ChargeScale = 0;
    #endif
    unsigned long long num  = 1000; //!< Arb, Number of particles to be collected before saving
    double TimeStepFactor   = 0.0005; //!< Arb, Multiplicative factor used to determine size of the timestep
    unsigned int Saves(100);    //!<Arb, Number of particles to be collected before saving in a run
    unsigned int reflectionsmax(50); //!< Arb, Number of reflections before rejecting particles


    // ************************************************** //


    // ***** RANDOM NUMBER GENERATOR        ***** //
    //!< Arb, Seed for the random number generator
    double seed     = 1.0;

    // ************************************************** //

    // ***** DETERMINE USER INPUT ***** //
    std::vector <std::string> sources;
    std::stringstream ss0;
    for (int i = 1; i < argc; ++i){ //!< Read command line input
        std::string arg = argv[i];
        int InputStatus(0);
        if     ( arg == "--help"    || arg == "-h" ){
            show_usage( argv[0]); return 0;
        }
        else if( arg == "--radius"  || arg == "-r" )
            InputStatus = InputFunction(argc,argv,i,ss0,Radius);
        else if( arg == "--spin"    || arg == "-s" )
            InputStatus = InputFunction(argc,argv,i,ss0,Spin);
        else if( arg == "--semix"   || arg == "-a1")
            InputStatus = InputFunction(argc,argv,i,ss0,a1);
        else if( arg == "--semiy"   || arg == "-a2")
            InputStatus = InputFunction(argc,argv,i,ss0,a2);
        else if( arg == "--semiz"   || arg == "-a3")
            InputStatus = InputFunction(argc,argv,i,ss0,a3);
        else if( arg == "--density"     || arg == "-d" )
            InputStatus = InputFunction(argc,argv,i,ss0,Density);
        else if( arg == "--efieldx"   || arg == "-ex" )
            InputStatus = InputFunction(argc,argv,i,ss0,EFieldx);
        else if( arg == "--efieldy"   || arg == "-ey" )
            InputStatus = InputFunction(argc,argv,i,ss0,EFieldy);
        else if( arg == "--efieldz"   || arg == "-ez" )
            InputStatus = InputFunction(argc,argv,i,ss0,EFieldz);
        else if( arg == "--potential"   || arg == "-p" )
            InputStatus = InputFunction(argc,argv,i,ss0,Potential);
        else if( arg == "--magfield"    || arg == "-m" )
            InputStatus = InputFunction(argc,argv,i,ss0,BMagIn);
        else if( arg == "--normalised"  || arg == "-n" )
            InputStatus = InputFunction(argc,argv,i,ss0,NormVars);
        else if( arg == "--etemp"   || arg == "-te")
            InputStatus = InputFunction(argc,argv,i,ss0,eTemp);
        else if( arg == "--edensity"    || arg == "-ne")
            InputStatus = InputFunction(argc,argv,i,ss0,eDensity);
        else if( arg == "--imass"   || arg == "-mi")
            InputStatus = InputFunction(argc,argv,i,ss0,iMass);
        else if( arg == "--itemp"   || arg == "-ti")
            InputStatus = InputFunction(argc,argv,i,ss0,iTemp);
        else if( arg == "--idensity"    || arg == "-ni")
            InputStatus = InputFunction(argc,argv,i,ss0,iDensity);
        else if( arg == "--ichance"     || arg == "-c" )
            InputStatus = InputFunction(argc,argv,i,ss0,iChance);
        else if( arg == "--zmaxcoeff"   || arg == "-u" )
            InputStatus = InputFunction(argc,argv,i,ss0,zMaxCoeff);
        else if( arg == "--zmincoeff"   || arg == "-l" )
            InputStatus = InputFunction(argc,argv,i,ss0,zMinCoeff);
        else if( arg == "--zboundforce" || arg == "-z" )
            InputStatus = InputFunction(argc,argv,i,ss0,ZBoundForce);
        else if( arg == "--impactpar"   || arg == "-b" )
            InputStatus = InputFunction(argc,argv,i,ss0,ImpactPar);
        else if( arg == "--forceimppar" || arg == "-f" )
            InputStatus = InputFunction(argc,argv,i,ss0,ForceImpPar);
        else if( arg == "--imax"    || arg == "-i" )
            InputStatus = InputFunction(argc,argv,i,ss0,imax);
        #if defined VARIABLE_CSCALE || defined VARIABLE_ASCALE
        else if( arg == "--jmin"    || arg == "-jm")
            InputStatus = InputFunction(argc,argv,i,ss0,jmin);
        else if( arg == "--jfin"    || arg == "-jf")
            InputStatus = InputFunction(argc,argv,i,ss0,jfin);
        #endif
        #if defined VARIABLE_CSCALE
        else if( arg == "--chargescale" || arg == "-cs")
            InputStatus = InputFunction(argc,argv,i,ss0,ChargeScale);
        #endif
        else if( arg == "--jmax"    || arg == "-j" )
            InputStatus = InputFunction(argc,argv,i,ss0,jmax);
        else if( arg == "--rmax"    || arg == "-rm")
            InputStatus = InputFunction(argc,argv,i,ss0,reflectionsmax);
        else if( arg == "--time"    || arg == "-t" )
            InputStatus = InputFunction(argc,argv,i,ss0,TimeStepFactor);
        else if( arg == "--number"  || arg == "-no")
            InputStatus = InputFunction(argc,argv,i,ss0,num);
        else if( arg == "--driftvel"    || arg == "-v" )
            InputStatus = InputFunction(argc,argv,i,ss0,DriftVel);
        else if( arg == "--seed"    || arg == "-se")
            InputStatus = InputFunction(argc,argv,i,ss0,seed);
        else if( arg == "--saves"   || arg == "-sa")
            InputStatus = InputFunction(argc,argv,i,ss0,Saves);
        else if( arg == "--output"  || arg == "-o" )
            InputStatus = InputFunction(argc,argv,i,ss0,suffix);
        #if defined CUSTOM_POTENTIAL
	else if( arg == "--customfield" || arg == "-cf")
	    InputStatus = InputFunction(argc,argv,i,ss0,input_custom_filename);
        #endif
        else{
            sources.push_back(argv[i]);
        }
        if( InputStatus == 1 ){
            std::cerr << "\nError! Input not recognised";
            return 1;
        }
    }
    // ************************************************** //

    // ***** INPUT CHECKING ***** //
    if( Radius < 0.0 )
        std::cerr << "\nError! Probe Radius is negative\nRadius : "
            << Radius;
    if( Density < 0.0 )
        std::cerr << "\nError! Probe Density is negative\nDensity : "
            << Density;
    if( eTemp < 0.0 )
        std::cerr << "\nError! Electron Temperature is negative\neTemp : "
            << eTemp;
    if( eDensity < 0.0 )
        std::cerr << "\nError! Electron Density is negative\neDensity : " << eDensity;
    if( iMass < 1.0 )
        std::cerr << "\nError! Ion Mass is less than 1\niMass : " << iMass;
    if( iTemp < 0.0 )
        std::cerr << "\nError! Ion Temperature is negative\niTemp : "
            << iTemp;
    if( iDensity < 0.0 )
        std::cerr << "\nError! Ion Density is negative\niDensity : "
            << iDensity;
    if( iChance < 0.0 && iChance != -0.5 )
        std::cout << "\nWarning! Chance of generating Ion is negative, "
            << "self-consistent flux assumed\niChance : " << iChance;
    if( zMinCoeff < 0.0 )
        std::cerr << "\nError! Lower Vertical Boundary Parameter is "
            << "negative\nzMinCoeff : " << zMinCoeff;
    if( zMaxCoeff < 0.0 )
        std::cerr << "\nError! Upper Vertical Boundary Parameter is "
            << "negative\nzMaxCoeff : " << zMaxCoeff;
    if( ZBoundForce < 0.0 )
        std::cerr << "\nError! Force Vertical Boundaries Parameter is "
            << "negative\nZBoundForce : " << ZBoundForce;
    if( ImpactPar < 0.0 )
        std::cerr << "\nError! Impact Parameter is negative\nImpactPar : "
            << ImpactPar;
    if( ForceImpPar < 0.0 )
        std::cerr << "\nError! Force Impact Parameter is negative\n"
            << "ForceImpPar : " << ForceImpPar;
    if( reflectionsmax < 0.0 )
        std::cerr << "\nError! Maximum Reflections Parameter is negative\n"
            << "reflectionsmax : " << reflectionsmax;
    if( imax < jmax )
        std::cerr << "\nError! Total particle goal less than captured particle "
            << "goal\nimax < jmax : "
            << imax << " < " << jmax;
    if( jmax < num )
        std::cout << "\nWarning! Save interval less than captured particle "
            << "goal. No Angular data recorded\nnum < jmax : "
            << num << " < " << jmax;
    #if defined VARIABLE_CSCALE || defined VARIABLE_ASCALE
    if( jfin == jmin ){
        std::cout << "\nWarning! Final collected equals minimum collected!\n"
            << "jfin == jmin : " << jfin << " == " << jmin << "\n"
            << "setting jfin = jmin+1";
        jfin = jmin+1;
    }
    #endif
    if( Saves > imax ){
        std::cout << "\nwarning! Saves greater than number of simulated "
            << "particles. Code won't run!\nsetting Saves = imax";
        Saves = imax;
    }
    #if defined CUSTOM_POTENTIAL
	std::string custom_filestring = FIELD_DIRECTORY+DIRECTORY_DELIM+input_custom_filename;
    // Check if file exists
    std::ifstream f(custom_filestring);
    if (!f.good()){
        std::cerr<<"\nError! Field file not found in "<<FIELD_DIRECTORY<<" directory.\n"
            << "Target File: "<<input_custom_filename<<"\n"
            << "Using Default File: "<<DEFAULT_FIELD_FILENAME<<"\n";
        input_custom_filename = DEFAULT_FIELD_FILENAME;
        //Note to self: some method for quitting the program. Put default string as constant defined elsewhere
    }
    #endif


    // ************************************************** //

    // ***** WRITE INPUT DATA FILE        ***** //
    //!< Save program configuration to a file with Input
    InputDataFile.open(filename + "_Input" + suffix);
    InputDataFile << "## Input Data File ##\n";
    InputDataFile << "#Input:\nr\ts\ta1\ta2\ta3\td\tp\tm\tn\tte\tne\tmi\tti\tni"
        << "\tc\tu\tl\tz\tb\tf\ti\tj\tt\tno\tv\tse\tsa\to";
    #if defined VARIABLE_CSCALE || defined VARIABLE_ASCALE
    InputDataFile << "\tjmin\tjfin";
    #endif
    #if defined CUSTOM_POTENTIAL
    InputDataFile << "\tcf";
    #endif
    InputDataFile << "\n" << Radius << "\t" << Spin << "\t" << a1 << "\t" << a2
        << "\t" << a3 << "\t" << Density  << "\t" << Potential << "\t" << BMagIn
        << "\t" << NormVars << "\t" << eTemp << "\t" << eDensity << "\t"
        << iMass << "\t" << iTemp << "\t" << iDensity << "\t" << iChance
        << "\t" << zMaxCoeff << "\t" << zMinCoeff << "\t" << ZBoundForce
        << "\t" << ImpactPar << "\t" << ForceImpPar << "\t" << imax << "\t"
        << jmax  << "\t" << TimeStepFactor << "\t" << num << "\t" << DriftVel
        << "\t" << seed << "\t" << Saves << "\t" << suffix;
    #if defined VARIABLE_CSCALE || defined VARIABLE_ASCALE
    InputDataFile << "\t" << jmin << "\t" << jfin;
    #endif
    #if defined CUSTOM_POTENTIAL
    InputDataFile << "\t" << input_custom_filename;
    #endif
    InputDataFile.close();

    #ifdef SAVE_TRACKS
    if( imax >= 1000 ){
        std::cout << "\nWarning! Attempting to store more than 1000 tracks!\n"
            << "i = " << imax;
        std::cout << "\nContinue?...\n";
        std::cin.get();
    }
    #endif

    // ************************************************** //

    // ***** ENACT NORMALISATION SCHEME ***** //

    //!< If species is positively charged, we assume it's a singly charged ion.
    //!< Otherwise, singly charged electron

    //!< If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
    double MASS = iMass*dimplconsts::Mp;  //!< kg, This is the Mass to which quantities are normalised
    a1 = 1.0/(a1*a1);
    a2 = 1.0/(a2*a2);
    a3 = 1.0/(a3*a3);
    double MassRatio    = sqrt(MASS/Me);
    double DustMass     = (4.0/3.0)*PI*pow(Radius,3)*Density;
    // If we're using S&L normalised units but iChance is undertermined
    if( NormVars ){
        if( iChance == 0.0 ){ // If we are simulating only Electrons
            // BMag normalised to Electrons
            BMag = sqrt(PI/2.0)*BMagIn*sqrt(Me*eTemp/echarge)/Radius;
        }else{  // If we are simulating only Ions or Ions and electrons.
            // BMag normalised to Ions
            BMag = sqrt(PI/2.0)*BMagIn*sqrt(MASS*iTemp/echarge)/Radius;
        }
        DriftVel = DriftVel*sqrt(echarge*iTemp/MASS);
    }else{
        // Convert from SI Potential to normalised potential
        BMag = BMagIn;
        Potential = Potential*echarge/(echarge*eTemp);
    }

    // Normalise TIME to the Gyro-Radius of an Ion at B=100T
    // Normalise MASS to Ion Mass
    // Normalise DISTANCE to Dust Radius
    // Normalise CHARGE to fundamental charge
    double MAGNETIC(1);
    double Tau = MASS/(echarge*MAGNETIC);

    // Normalised Charge,
    double PotentialNorm
        = Potential*(eTemp*echarge)*4*PI*epsilon0*Radius/pow(echarge,2.0);
    double DriftNorm
        = DriftVel*Tau/(Radius);
    double DebyeLength
        = sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2.0)))/Radius;
    double A_Coulomb
        = MASS/(4.0*PI*epsilon0*MAGNETIC*MAGNETIC*Radius*Radius*Radius);

    #ifdef VARIABLE_CSCALE
    if( ChargeScale == 0.0 ){
	ChargeScale = eTemp*4.0*PI*epsilon0*Radius/(2.0*echarge); // Divide by 2 as we don't want full scale
    }
    #endif
    #ifdef VARIABLE_ASCALE
    double AngularScalei
        = MASS*iDensity*sqrt(echarge*iTemp/MASS)*Tau*0.001/(Radius);
    double AngularScalee
        = MASS*eDensity*sqrt(echarge*eTemp/MASS)*Tau*0.001/(MassRatio*Radius);
    #else
    double AngularScalei = 1.0;
    #endif

    // ************************************************** //


    // ***** DEFINE FIELD PARAMETERS        ***** //

    threevector EField_Background(EFieldx,EFieldy,EFieldz);
    threevector Bhat(0.0,0.0,1.0);  // Direction of magnetic field, z dir.

    // ************************************************** //

    // ***** TRIGGER CREATION OF CUSTOM ELECTRIC FIELD MAP STORAGE ***** //
    #ifdef CUSTOM_POTENTIAL
	custom_filestring = FIELD_DIRECTORY + DIRECTORY_DELIM + input_custom_filename;
    bool is_spherical_injection = false;
    bool is_cylindrical_injection = false;
    #ifdef SPHERICAL_INJECTION
    is_spherical_injection = true;
    #else
    is_cylindrical_injection = true;
    #endif
    Field_Map this_field_map(custom_filestring, Radius, DebyeLength, is_spherical_injection, is_cylindrical_injection, a1, a2, a3);
    
    #endif
    // ************************************************** //


    // ***** DEFINE SIMULATION SPACE        ***** //
    // See: https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    // For reasoning on choice of Velocity vector
    // Normalised Ion/Electron Thermal velocity
    double iThermalVel  = sqrt(echarge*iTemp/MASS)*(Tau/Radius);
    double eThermalVel  = sqrt(echarge*eTemp/Me)*(Tau/Radius);

    double iRhoTherm = 0.0;  // Ion gyro-radii are zero by Default
    double eRhoTherm = 0.0;  // Electron gyro-radii are zero by Default

    double TimeStepe = TimeStepFactor;
    double TimeStepi = TimeStepFactor;
    if( BMag == 0.0 && Potential == 0.0 ){
        TimeStepi = 0.01/iThermalVel;
        TimeStepe = 0.01/eThermalVel;
    }
    if( BMag != 0.0 ){ // If Magnetic field is non-zero
        // Calculate thermal GyroRadius for ions and electrons normalised to
        // dust grain radii
        iRhoTherm   = iThermalVel/(BMag/MAGNETIC);
        eRhoTherm   = eThermalVel/(pow(MassRatio,2)*BMag/MAGNETIC);

        if( NormVars ){
            iRhoTherm   = 1.0/BMagIn;
            eRhoTherm   = 1.0/BMagIn;

            // If there is a finite probability of simulating ions
            if( iChance != 0.0 ){
                eRhoTherm = 1.0/(BMagIn*MassRatio);
            }
        }
        // Time step limited by gyro-motion
        if( TimeStepFactor*eRhoTherm/eThermalVel < TimeStepe ){
            // 1% of Gyro-radius size
            TimeStepe = TimeStepFactor*eRhoTherm/eThermalVel;
        }
        if( TimeStepFactor*iRhoTherm/iThermalVel < TimeStepi ){
            // 1% of Gyro-radius size
            TimeStepi = TimeStepFactor*iRhoTherm/iThermalVel;
        }
        if( Potential == 0.0 ){
            // Uncharged sphere, time step limited by gyro-motion or thermal
            // velocity
            // 1% of Gyro-radius size
            TimeStepe = TimeStepFactor*eRhoTherm/eThermalVel;
            TimeStepi = TimeStepFactor*iRhoTherm/iThermalVel;
        }
        // Check timestep less than dust radius
        if( TimeStepFactor/iThermalVel < TimeStepi ){
            TimeStepi   = TimeStepFactor/iThermalVel;
        }
        if( TimeStepFactor/eThermalVel < TimeStepe ){
            TimeStepe   = TimeStepFactor/eThermalVel;
        }
	// Check timestep less than sheath size
        #ifdef DEBYE_POTENTIAL
        if( TimeStepFactor*DebyeLength/iThermalVel < TimeStepi ){
            TimeStepi   = TimeStepFactor*DebyeLength/iThermalVel;
        }
        if( TimeStepFactor*DebyeLength/eThermalVel < TimeStepe ){
            TimeStepe   = TimeStepFactor*DebyeLength/eThermalVel;
        }
        #endif
        #ifdef PARABOLIC_POTENTIAL
        if( TimeStepFactor*DebyeLength/iThermalVel < TimeStepi ){
            TimeStepi   = TimeStepFactor*DebyeLength/iThermalVel;
        }
        if( TimeStepFactor*DebyeLength/eThermalVel < TimeStepe ){
            TimeStepe   = TimeStepFactor*DebyeLength/eThermalVel;
        }
        #endif
    }

    // Account for non-spherical impact factor, chose larger of two semi axis
    double semiaxisDistortion=1.0;
    if( a1 < a2 ){
        semiaxisDistortion=a1;
    }else{
        semiaxisDistortion=a2;
    }
    double iImpactParameter
        = 1.0/sqrt(semiaxisDistortion)+ImpactPar*(iRhoTherm+DebyeLength);
    double eImpactParameter
        = 1.0/sqrt(semiaxisDistortion)+ImpactPar*(eRhoTherm+DebyeLength);

    #ifdef SELF_CONS_CHARGE
    #else
    if( PotentialNorm == 0.0 ){
        iImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*iRhoTherm;
        eImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*eRhoTherm;
    }
    #endif
    
    #ifdef CUSTOM_POTENTIAL
    const double field_map_rho_max = this_field_map.get_map_rho_max();
    // If -f not used, scale as below
    iImpactParameter = 3*(iRhoTherm+DebyeLength);
    eImpactParameter = 3*(eRhoTherm);
    if (iImpactParameter>field_map_rho_max||eImpactParameter>field_map_rho_max){
        std::cerr<<"ERROR: the custom field map domain is not large enough in the radial direction.\n\tPlease choose a larger coverage of the custom field map."<<std::endl;
        //Note to self: quit program
    }
    if( ForceImpPar > 0.0 ){
        if (ForceImpPar > field_map_rho_max){
            std::cerr<<"ERROR: zboundforce provided is larger than the custom field map domain.\n\tPlease redefine the zboundforce or choose a larger coverage of the custom field map."<<std::endl;
            //Note to self: quit program
        } else {
            iImpactParameter = ForceImpPar;
            eImpactParameter = ForceImpPar;
            //Note to self: shrink map accordingly
        }
    }
    #else
    if( ForceImpPar > 0.0 ){
        iImpactParameter = ForceImpPar;
        eImpactParameter = ForceImpPar;
    }
    #endif
    

    double ezmax = 1.0/sqrt(a3);
    double ezmin = -1.0/sqrt(a3);
    double izmax = 1.0/sqrt(a3);
    double izmin = -1.0/sqrt(a3);

    #ifdef COULOMB_POTENTIAL
    // Balance Coulomb to thermal energy
    double iCoulombImpactParameter
        = 10.0*fabs(echarge*echarge*PotentialNorm
        /(2*PI*epsilon0*echarge*iTemp))/Radius;
    double eCoulombImpactParameter
        = 10.0*fabs(echarge*echarge*PotentialNorm
        /(2*PI*epsilon0*echarge*eTemp))/Radius;
    ezmax += zMaxCoeff*eCoulombImpactParameter;
    ezmin -= zMinCoeff*eCoulombImpactParameter;
    izmax += zMaxCoeff*iCoulombImpactParameter;
    izmin -= zMinCoeff*iCoulombImpactParameter;
    #elif defined DEBYE_POTENTIAL
    ezmax += zMaxCoeff*DebyeLength;
    ezmin -= zMinCoeff*DebyeLength;
    izmax += zMaxCoeff*DebyeLength;
    izmin -= zMinCoeff*DebyeLength;
    #elif defined PARABOLIC_POTENTIAL
    ezmax += zMaxCoeff*DebyeLength;
    ezmin -= zMinCoeff*DebyeLength;
    izmax += zMaxCoeff*DebyeLength;
    izmin -= zMinCoeff*DebyeLength;
    #elif defined CUSTOM_POTENTIAL
    ezmax = this_field_map.get_map_z_max();
    ezmin = this_field_map.get_map_z_min();
    izmax = ezmax;
    izmin = ezmin;
    #else
    ezmax += zMaxCoeff;
    ezmin -= zMinCoeff;
    izmax += zMaxCoeff;
    izmin -= zMinCoeff;
    #endif
    
    #ifdef CUSTOM_POTENTIAL
    if(ZBoundForce>ezmax){
        std::cerr<<"ERROR: zboundforce provided is larger than the custom field map domain.\n\tPlease redefine the zboundforce or choose a larger coverage of the custom field map."<<std::endl;
        //Note to self: quit program
    } else {
        ezmax = ZBoundForce;
        ezmin = -ezmax;
        izmax = ezmax;
        izmin = ezmin;
        //Note to self: shrink map accordingly.
    }
    #else
    if( ZBoundForce > 0.0 ){
        ezmax = 1.0/sqrt(a3)+ZBoundForce;
        ezmin = -1.0/sqrt(a3)-ZBoundForce;
        izmax = ezmax;
        izmin = ezmin;
    }
    #endif

    // ************************************************** //

    // ***** SHRINK CUSTOM ELECTRIC FIELD MAP ACCORDINGLY ***** //
    #ifdef CUSTOM_POTENTIAL
    this_field_map.shrink_map_to_fit(std::min(ezmin, izmin), std::max(ezmax, izmax), std::max(eImpactParameter, iImpactParameter));
    #endif
    // ************************************************** //

    // ***** CONFIGURE SPHERICAL INJECTION        ***** //
    // find the maximum value of the theta pdf
    #ifdef SPHERICAL_INJECTION
    double iThetaPDFMax = thetaPDFMax(DriftNorm,iThermalVel);
    double eThetaPDFMax = thetaPDFMax(DriftNorm,eThermalVel);
    #else
    double iThetaPDFMax = 0.0;
    double eThetaPDFMax = 0.0;
    #endif

    // ************************************************** //

    // ***** DEFINE PROBABILITY OF ION GENERATION   ***** //
    // Define ratio of flux of electrons to ions
    // The min and max heights must match as long as the ProbabilityOfIon
    // is the same for both the top and bottom surfaces.
    assert(fabs(ezmax)==fabs(ezmin));
    assert(fabs(izmax)==fabs(izmin));
    #ifdef BOLTZMANN_DENSITY
        #ifdef SPHERICAL_INJECTION
            #ifdef COULOMB_POTENTIAL
            double BoltzmanneDensity
                = eDensity*exp(PotentialNorm*echarge*echarge
                /(4.0*PI*epsilon0*eImpactParameter*Radius*echarge*eTemp));
            double BoltzmanniDensity
                = iDensity*exp(-PotentialNorm*echarge*echarge
                /(4.0*PI*epsilon0*iImpactParameter*Radius*echarge*iTemp));
            #else
            double BoltzmanneDensity
                = eDensity*exp(PotentialNorm*echarge*echarge
                *exp(-(eImpactParameter-1.0)/DebyeLength)
                /(4.0*PI*epsilon0*eImpactParameter*Radius*echarge*eTemp));
            double BoltzmanniDensity
                = iDensity*exp(-PotentialNorm*echarge*
                *exp(-(iImpactParameter-1.0)/DebyeLength)
                /(4.0*PI*epsilon0*iImpactParameter*Radius*echarge*iTemp));
            #endif
        #else
            #ifdef COULOMB_POTENTIAL
            double BoltzmanneDensity
                = eDensity*exp(PotentialNorm*echarge*echarge
                /(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
            double BoltzmanniDensity
                = iDensity*exp(-PotentialNorm*echarge*echarge
                /(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
            #else
            double BoltzmanneDensity
                = eDensity*exp(PotentialNorm*echarge*echarge
                *exp(-(ezmax-1.0)/DebyeLength)
                /(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
            double BoltzmanniDensity
                = iDensity*exp(-PotentialNorm*echarge*echarge
                *exp(-(izmax-1.0)/DebyeLength)
                /(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
            #endif
        #endif
    #else
        double BoltzmanneDensity = eDensity;
        double BoltzmanniDensity = iDensity;
    #endif

    double PosFluxi
        = BoltzmanniDensity*((iThermalVel/sqrt(2.0*PI))
        *exp(-0.5*DriftNorm*DriftNorm/(iThermalVel*iThermalVel))
        +DriftNorm*0.5*(1.0+erf(DriftNorm/(sqrt(2.0)*iThermalVel))))
        *pow(iImpactParameter,2);
    double NegFluxi
        = BoltzmanniDensity*((iThermalVel/sqrt(2.0*PI))
        *exp(-0.5*DriftNorm*DriftNorm/(iThermalVel*iThermalVel))
        -DriftNorm*0.5*(1.0+erf(-DriftNorm/(sqrt(2.0)*iThermalVel))))
        *pow(iImpactParameter,2);
    double PosFluxe
        = BoltzmanneDensity*((eThermalVel/sqrt(2.0*PI))
        *exp(-0.5*DriftNorm*DriftNorm/(eThermalVel*eThermalVel))
        +DriftNorm*0.5*(1.0+erf(DriftNorm/(sqrt(2.0)*eThermalVel))))
        *pow(eImpactParameter,2);
    double NegFluxe
        = BoltzmanneDensity*((eThermalVel/sqrt(2.0*PI))
        *exp(-0.5*DriftNorm*DriftNorm/(eThermalVel*eThermalVel))
        -DriftNorm*0.5*(1.0+erf(-DriftNorm/(sqrt(2.0)*eThermalVel))))
        *pow(eImpactParameter,2);

    double TotalFluxProbability = PosFluxi + NegFluxi + PosFluxe + NegFluxe;

    double ProbOfPosFluxe = PosFluxe / TotalFluxProbability;
    double ProbOfNegFluxe = NegFluxe / TotalFluxProbability;
    double ProbOfPosFluxi = PosFluxi / TotalFluxProbability;
    double ProbOfNegFluxi = NegFluxi / TotalFluxProbability;

    double ProbabilityOfIon = ProbOfPosFluxi+ProbOfNegFluxi;
    if( iChance >= 0.0 && iChance <= 1.0 ){
        ProbabilityOfIon = iChance;
        TotalFluxProbability
            = ProbabilityOfIon*(PosFluxi + NegFluxi)
            + (1.0-ProbabilityOfIon)*(PosFluxe + NegFluxe);

        ProbOfPosFluxe = (1.0-ProbabilityOfIon)*PosFluxe / TotalFluxProbability;
        ProbOfNegFluxe = (1.0-ProbabilityOfIon)*NegFluxe / TotalFluxProbability;
        ProbOfPosFluxi = ProbabilityOfIon*PosFluxi / TotalFluxProbability;
        ProbOfNegFluxi = ProbabilityOfIon*NegFluxi / TotalFluxProbability;
    }

    // ************************************************** //


    // ***** SEED RANDOM NUMBER GENERATOR IN THREADS***** //
    std::random_device rd;      // Create Random Device
    std::vector<std::mt19937> randnumbers;
    std::uniform_real_distribution<double> rad(0, 1);
    for(int p = 0; p < omp_get_max_threads(); p ++){
        randnumbers.push_back(std::mt19937(seed+p));
    }

    // ************************************************** //


    // ***** OPEN DATA FILES WITH HEADER         ***** //
    time_t now = time(0); // Get the time of simulation
    char * dt = ctime(&now);

    OPEN_CHARGE();  HEAD_CHARGE();
    OPEN_MOM();     HEAD_MOM();
    OPEN_AVEL();    HEAD_AVEL();
    OPEN_LMOM();    HEAD_LMOM();
    OPEN_CHA();     HEAD_CHA();
    OPEN_SPOS();    HEAD_SPOS();
    OPEN_EPOS();    HEAD_EPOS();
    OPEN_APP();     HEAD_APP()
    OPEN_CURR();    HEAD_CURR();
    OPEN_TOT();     HEAD_TOT();
    OPEN_REF();     HEAD_REF();
    OPEN_AXIS();    HEAD_AXIS();

    unsigned long long NumberOfIons = ProbabilityOfIon*imax;
    unsigned long long NumberOfElectrons = imax-NumberOfIons;

    RunDataFile.open(filename + suffix);
    RunDataFile << "## Run Data File ##\n";
    RunDataFile << "#Date: " << dt;
    RunDataFile << "#Input:\t\tValue\n\nimax (arb #):\t\t"<<imax
        <<"\njmax (arb #):\t\t"<<jmax<<"\nIon # (arb #):\t\t" << NumberOfIons
        <<"\nElec # (arb #):\t\t"<<NumberOfElectrons<<"\nProbOfIon (arb):\t"
        <<ProbabilityOfIon<<"\n\nElectron Gyro (1/Radius):\t"<<eRhoTherm
        <<"\nElectron Temp (eV):\t\t"<<eTemp<<"\nElec Density (m^-^3):\t\t"
        <<BoltzmanneDensity<<"\nElectron IP (1/Radius):\t\t"<<eImpactParameter
        <<"\nElectron zmax (1/Radius):\t"<<ezmax
        <<"\nElectron zmin (1/Radius):\t"<<ezmin<<"\nElecs Timestep (Tau):\t\t"
        <<TimeStepe<<"\nProbOfPosFluxe (arb):\t\t"<<ProbOfPosFluxe
        <<"\nProbOfNegFluxe (arb):\t\t"<<ProbOfNegFluxe
        <<"\n\nIon Gyro (1/Radius):\t"<<iRhoTherm<<"\nIon Temp (eV):\t\t"
        <<iTemp<<"\nIon Density (m^-^3):\t"<<BoltzmanniDensity
        <<"\nIon IP (1/Radius):\t"<<iImpactParameter<<"\nIon zmax (1/Radius):\t"
        <<izmax<<"\nIon zmin (1/Radius):\t"<<izmin<<"\nIon Timestep (Tau):\t"
        <<TimeStepi<<"\nProbOfPosFluxi (arb):\t"<<ProbOfPosFluxi
        <<"\nProbOfNegFluxi (arb):\t"<<ProbOfNegFluxi<<"\n\nRadius (m):\t\t"
        <<Radius<<"\nSpin (1/Tau):\t\t"<<Spin*Tau<<"\na1 (1/Radius):\t\t"
        <<(1.0/sqrt(a1))<<"\na2 (1/Radius):\t\t"<<(1.0/sqrt(a2))
        <<"\na3 (1/Radius):\t\t"<<(1.0/sqrt(a3))<<"\nDensity (kg m^-^3):\t"
        <<Density<<"\nCharge (1/echarge):\t\t"<<PotentialNorm
        <<"\nB Field (T or Radius/GyroRad):\t"<<BMag
        <<"\nDebyeLength (1/Radius):\t\t"<<DebyeLength
        <<"\nDrift Norm (Radius/Tau):\t"<<DriftNorm
        <<"\nTime Norm [Tau] (s):\t\t"<<Tau<<"\n\n"<<"RNG Seed (arb):\t\t"
        <<seed<<"\nOMP_THREADS (arb):\t"<<omp_get_max_threads()<<"\n\n";

    #ifdef SPHERICAL_INJECTION
        RunDataFile << "* SPHERICAL INJECTION *\n";
    #else
        RunDataFile << "* CYLINDRICAL INJECTION *\n";
    #endif
    #ifdef NO_SPHERE
        RunDataFile << "* NO SPHERE *\n";
    #elif defined DISK
        RunDataFile << "* DISK *\n";
    #endif

    #ifdef DEBYE_POTENTIAL
        RunDataFile << "* DEBYE HUCKEL POTENTIAL *\n\n";
    #endif
    #ifdef PARABOLIC_POTENTIAL
        RunDataFile << "* PARABOLIC POTENTIAL *\n\n";
    #endif

    #ifdef COULOMB_POTENTIAL
        RunDataFile << "* COULOMB POTENTIAL *\n\n";
    #endif
    #ifdef CUSTOM_POTENTIAL
	RunDataFile << "* CUSTOM POTENTIAL *\n\n";
    #endif
    #ifdef BOLTZMANN_DENSITY
        RunDataFile << "* BOLTZMANN DENSITY *\n\n";
    #endif

    // ************************************************** //

    // ***** DEFINE SHARED VARIABLES USED GLOBALLY IN OPERATION    ***** //
    threevector TotalAngularVel(0.0,0.0,Spin*Tau);
    #ifdef FREE_AXIS
        threevector AxisAVel(0.0,0.0,Spin*Tau); //!< Dust axis velocity
    #endif

    threevector TotalAngularMom(0.0,0.0,0.0);
    threevector TotalInjectedMomCaptured(0.0,0.0,0.0);
    threevector TotalInjectedMomMissed(0.0,0.0,0.0);
    threevector TotalLostMom(0.0,0.0,0.0);

    #ifdef VARIABLE_CSCALE
    // Initialise the mean charge as starting charge
    double MeanChargeSave    = PotentialNorm;
    double TotalChargeInSave = PotentialNorm;
    double MeanChargeDiff(0.0);    // Initial Mean charge diff is zero
    double OldMeanChargeDiff(0.0); // Initial Mean charge diff is zero
    #endif

    #ifdef VARIABLE_ASCALE
    threevector MeanAngularVelSave(0.0,0.0,Spin*Tau);
    threevector TotalAngularVelThisStep(0.0,0.0,Spin*Tau);
    threevector MeanAngularVelDiff(0.0,0.0,0.0);
    threevector OldMeanAngularVelDiff(0.0,0.0,0.0);
    #endif

    #ifdef SAVE_MOM
    threevector LinearMomentumSum(0.0,0.0,0.0);
    threevector AngularMomentumSum(0.0,0.0,0.0);
    #endif

    #ifdef TEST_ANGMOM
    threevector INITIAL_AMOM(0.0,0.0,0.0);
    threevector FINAL_AMOM(0.0,0.0,0.0);
    #endif

    unsigned long long j(0), i(0), RegeneratedParticles(0), TrappedParticles(0);
    unsigned long long MissedParticles(0), TotalNum(0);
    unsigned long long j_ThisSave(0), e_simulated(0), i_simulated(0);
    long long CapturedCharge(0), RegeneratedCharge(0), TrappedCharge(0);
    long long MissedCharge(0), TotalCharge(0);

    unsigned long long smax = imax / Saves;

    // ************************************************** //


    // ***** BEGIN LOOP OVER NUMBER OF SAVES    ***** //

    for( unsigned int s=1; s <= smax; s ++){ // Loop over the number of saves
        unsigned long long IMAX = (s*imax)/smax;

        // ***** REOPEN DATA FILES          ***** //
        RunDataFile.close();
        RunDataFile.clear();
        RunDataFile.open(filename+suffix, std::fstream::app);

        REOPEN_CHARGE();
        REOPEN_MOM();
        REOPEN_AVEL();
        REOPEN_LMOM();
        REOPEN_CHA();
        REOPEN_SPOS();
        REOPEN_EPOS();
        REOPEN_APP();
        REOPEN_CURR();
        REOPEN_TOT();
        REOPEN_REF();
        REOPEN_AXIS();

        // ************************************************** //

        // ***** BEGIN PARALLELISATION REGION    ***** //
        #pragma omp parallel shared(TotalAngularVel,TotalAngularMom,TotalInjectedMomCaptured,TotalInjectedMomMissed,TotalLostMom,j) PRIVATE_FILES()
        {
        // ***** DEFINE PRIVATE VARIABLES USED IN THREAD OPERATION    ***** //
        threevector Position(0.0,0.0,0.0);
        threevector Velocity(0.0,0.0,0.0);
        threevector BField(0.0,0.0,0.0);
        threevector EField(0.0,0.0,0.0);
        threevector OldPosition(0.0,0.0,0.0);
        threevector AngularVel(0.0,0.0,0.0);
        threevector InitialPos(0.0,0.0,0.0);
        threevector InitialVel(0.0,0.0,0.0);
        threevector FinalPosition(0.0,0.0,0.0);
        threevector AngularMom(0.0,0.0,0.0);

        threevector DummyPosition(0.0,0.0,0.0);

        unsigned int reflections(0);

        double BMagNorm;
        double ImpactParameter;
        double ThermalVel;
        double xThetaPDFMax;
        double zmax;
        double zmin;
        double TimeStep;
        double TotalTime;
        double SpeciesMass;
        double ProbUpper;

        int SPEC_CHARGE=-1;

        bool EdgeCondition = true;
        bool SphereCondition = true;

        // ************************************************** //


        // ***** BEGIN PARALLELISED LOOP OVER PARTICLES    ***** //
        #pragma omp for
        // Loop over maximum number of particles to generate
        for( i=(IMAX-imax/smax); i < IMAX; i ++){
            // ***** LOOP UNTIL JMAX PARTICLES ARE COLLECTED    ***** //
            // Loop until we reach a certain number of particles jmax
            if( j <= jmax ){

                //std::cout << "\n" << omp_get_thread_num() << "/"
                //<< omp_get_num_threads();

                // ***** DETERMINE IF IT'S AN ELECTRON OR ION ***** //
                // MassRatio squared because of absence of mass in Boris solver
                #pragma omp critical
                {
                if( rad(randnumbers[omp_get_thread_num()]) < ProbabilityOfIon
                    && i_simulated < NumberOfIons ){
                    // If this is the case, we need to generate an ion
                    BMagNorm = BMag/MAGNETIC;
                    ImpactParameter=iImpactParameter;
                    ThermalVel=iThermalVel;
                    xThetaPDFMax=iThetaPDFMax;
                    zmax    = izmax;
                    zmin    = izmin ;
                    TimeStep = TimeStepi;
                    SpeciesMass = 1.0;
                    ProbUpper = ProbOfPosFluxi/(ProbOfPosFluxi+ProbOfNegFluxi);
                    SPEC_CHARGE=1;
                    i_simulated = i_simulated + 1;

                }else{ // If this is the case, we need to generate an electron
                    if( e_simulated < NumberOfElectrons ){
                        BMagNorm = BMag*pow(MassRatio,2)/MAGNETIC;
                        ImpactParameter=eImpactParameter;
                        ThermalVel=eThermalVel;
                        xThetaPDFMax=eThetaPDFMax;
                        zmax= ezmax; // Top of Simulation Domain, in Dust Radii
                        zmin= ezmin; // Top of Simulation Domain, in Dust Radii
                        TimeStep = TimeStepe;
                        SpeciesMass = 1.0/pow(MassRatio,2);
                        ProbUpper
                            = ProbOfPosFluxe/(ProbOfPosFluxe+ProbOfNegFluxe);
                        SPEC_CHARGE=-1;
                        e_simulated = e_simulated + 1;
                    }else if( i_simulated < NumberOfIons ){
                        BMagNorm = BMag/MAGNETIC;
                        ImpactParameter=iImpactParameter;
                        ThermalVel=iThermalVel;
                        xThetaPDFMax=iThetaPDFMax;
                        zmax    = izmax;
                        zmin    = izmin ;
                        TimeStep = TimeStepi;
                        SpeciesMass = 1.0;
                        ProbUpper
                            = ProbOfPosFluxi/(ProbOfPosFluxi+ProbOfNegFluxi);
                        SPEC_CHARGE=-1;
                        SPEC_CHARGE=1;
                        i_simulated = i_simulated + 1;
                    }else{
                        std::cerr << "\nWARNING! EXCESS PARTICLES SIMULATED";
                        zmax    = 1000.0;
                        zmin    = 999.0;
                        ImpactParameter = 1.0;
                    }
                }
                BField = BMagNorm*Bhat;
                }

                // ************************************************** //


                // ***** GENERATE AN ORBIT ***** //
                GenerateOrbit(Position,Velocity,ImpactParameter,ProbUpper,zmin,
                    zmax,DriftNorm,ThermalVel,xThetaPDFMax,
                    randnumbers[omp_get_thread_num()]);

                TotalTime = 0.0;
                InitialPos = Position;
                InitialVel = Velocity;

                // ************************************************** //


                // ***** TESTING AREA               ***** //
                // ***** ANGULAR-MOMENTUM TEST          ***** //
                #ifdef TEST_ANGMOM
                #pragma omp critical
                {
                    // For Angular Momentum Calculations
                    ADD_I_AMOM(SpeciesMass*(Position^Velocity));
                    PRINT_AMOM("IMom = "); PRINT_AMOM(INITIAL_AMOM);
                    PRINT_AMOM("\n");
                }
                #endif
                // ***** VELOCITY-POSITION DISTRIBUTION TEST    ***** //
                #ifdef TEST_VELPOSDIST
                #pragma omp critical
                {
                    // For debugging look at initial positions
                    PRINT_VPD(Position); PRINT_VPD("\t");
                    PRINT_VPD(Velocity); PRINT_VPD("\t");
                    PRINT_VPD(Velocity*(Position.getunit())); PRINT_VPD("\t");
                    PRINT_VPD( sqrt(pow(Velocity.getx(),2)
                        +pow(Velocity.gety(),2))*SpeciesMass);
                    PRINT_VPD("\n");    // For debugging look at gyro-radii
                }
                #endif
                // ***** ENERGY TEST: MEASURE INITIAL ENERGY    ***** //
                C_INITIAL_VEL(); D_INITIAL_VEL();   // For energy calculations
                C_INITIAL_POT(); D_INITIAL_POT();   // For energy calculations
                #ifdef SAVE_APPROACH
                    double MinPos
                        = sqrt(zmax*zmax+ImpactParameter*ImpactParameter);
                #endif
                // ***** RECORD TRACK DATA, DEBUG AND TEST  ***** //
                OPEN_TRACK(filename + "_Track_" + std::to_string(i) + suffix);
                RECORD_TRACK(Position);
                RECORD_TRACK("\t"); RECORD_TRACK(Velocity);
                RECORD_TRACK("\t"); RECORD_TRACK(sqrt(
                    Position.getx()*Position.getx()
                    +Position.gety()*Position.gety()
                    +Position.getz()*Position.getz()));

                // ************************************************** //


                // ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
                // Calculate Electric Field
                #ifdef DEBYE_POTENTIAL
                    EField = DebyeHuckelField(Position,PotentialNorm,
                        DebyeLength,A_Coulomb)+EField_Background;
                #endif
                #ifdef PARABOLIC_POTENTIAL
                    EField = ParabolicField(Position,PotentialNorm,
                        DebyeLength,A_Coulomb)+EField_Background;
                #endif
                #ifdef COULOMB_POTENTIAL
                    EField = CoulombField(Position,PotentialNorm,A_Coulomb)
                        +EField_Background;
                #endif
                #ifdef CUSTOM_POTENTIAL
                    EField = CustomField(Position,this_field_map)
                        +EField_Background;
                #endif

                UpdateVelocityBoris(SpeciesMass,EField,BField,-0.5*TimeStep,
                    Velocity,SPEC_CHARGE);

                OldPosition.setx(0.0);
                OldPosition.sety(0.0);
                OldPosition.setz(0.0);

                // ************************************************** //


                // ***** DEFINE SIMULATION BOUNDARY CONDITIONS       ***** //
                #ifdef SPHERICAL_INJECTION
                    EdgeCondition = ((Position.getx()*Position.getx()
                        +Position.gety()*Position.gety()
                        +Position.getz()*Position.getz())
                    <= ImpactParameter*ImpactParameter*1.01);
                #else
                    EdgeCondition = (Position.getz() >= zmin
                        && Position.getz() <= zmax);
                #endif


                #ifdef NO_SPHERE
                    SphereCondition = true;
                #elif defined DISK
                    SphereCondition = ( (Velocity.getz()*Position.getz() < 0.0)
                        && (sqrt(Position.getx()*Position.getx()*a1
                        +Position.gety()*Position.gety()*a2) < 1.0) );
                #else
                    SphereCondition = (sqrt(Position.getx()*Position.getx()*a1
                        +Position.gety()*Position.gety()*a2
                        +Position.getz()*Position.getz()*a3) > 1.0);
                #endif
                // ************************************************** //

                // ***** DEFINE REFLECTION TERMINATION CONDITION       ***** //
                // While we don't exceed a specified number of reflections to
                // catch trapped orbits AND while the particle is not inside
                // the sphere and not outside the simulation domain
                reflections=0;
                unsigned int r_max = reflectionsmax;
                // In this case, potential is repulsive
                if( PotentialNorm*SPEC_CHARGE > 0.0 ){
                    r_max = 1;
                    #ifdef MAGNETIC_BOTTLE
                    #ifdef COULOMB_POTENTIAL
                    double vx = Velocity.getx();
                    double vy = Velocity.gety();
                    double vz = Velocity.getz();

                    double r0 = sqrt(Position.getx()*Position.getx()
                        +Position.gety()*Position.gety());

                    double delta_a = 2.0*sqrt((vx*vx+vy*vy)+vz*vz/(2.0*PI)
                        +2*PotentialNorm/(BMag*BMag));
                    //std::cout << "\ne-\tr0 = " << r0 << "\t\tdelta_a = " << delta_a << "\n";
                    if(  delta_a < 0.0 || delta_a >= 2.0 ||
                        delta_a/2 < r0 ){
                        EdgeCondition = false;
                    }
                    #endif
                    #endif
/*                  // Compare vertical kinetic energy to potential.
                    // it's smaller, reject orbit
                    // immediately before simulating
                    #ifdef DEBYE_POTENTIAL
                    if( fabs((PotentialNorm*echarge*echarge
                        /(4.0*PI*epsilon0*Radius))
                        *(1.0-exp((zmax/DebyeLength)-(1.0/DebyeLength))/zmax))
                        > 0.5*SpeciesMass*Mp*Velocity.getz()*Velocity.getz()
                        *Radius*Radius/(Tau*Tau) ){
                            EdgeCondition = false;
                        }
                    #endif
                    #ifdef COULOMB_POTENTIAL
                        if( fabs((PotentialNorm*echarge*echarge
                            /(4.0*PI*epsilon0*Radius))*(1.0-(1.0/zmax)))
                            > (0.5*SpeciesMass*Mp*Velocity.getz()
                            *Velocity.getz()*Radius*Radius/(Tau*Tau)) ){
                            EdgeCondition = false;
                        }
                    #endif*/
                }else{
                    #ifdef MAGNETIC_BOTTLE
                    #ifdef COULOMB_POTENTIAL
                    double vx = Velocity.getx();
                    double vy = Velocity.gety();
                    double vz = Velocity.getz();
                    //double r0 = sqrt(Position.getx()*Position.getx()
                    //    +Position.gety()*Position.gety());

                    double delta_a = 2.0*sqrt((vx*vx+vy*vy)+vz*vz/(2.0*PI)
                        +2*PotentialNorm/(BMag*BMag));
                    //std::cout << "\ni+\tr0 = " << r0 << "\t\tdelta_a = " << delta_a << "\n";
                    if( delta_a <= 2.0 ){
                        EdgeCondition = true;
                    }else if( 1 <= sqrt(delta_a-1.0) ){
                        EdgeCondition = true;
                    }else{
                        EdgeCondition = false;
                    }
                    #endif
                    #endif
		}
                // ************************************************** //


                // ***** DO PARTICLE PATH INTEGRATION       ***** //
                while( SphereCondition && EdgeCondition
                    && reflections < r_max ){

                    #ifdef DEBYE_POTENTIAL
                        EField = DebyeHuckelField(Position,PotentialNorm,
                            DebyeLength,A_Coulomb)+EField_Background;
                    #endif
                    #ifdef PARABOLIC_POTENTIAL
                        EField = ParabolicField(Position,PotentialNorm,
                            DebyeLength,A_Coulomb)+EField_Background;
                    #endif
                    #ifdef COULOMB_POTENTIAL
                        EField = CoulombField(Position,PotentialNorm,
                            A_Coulomb)+EField_Background;
                    #endif
                    #ifdef CUSTOM_POTENTIAL
                        EField = CustomField(Position,this_field_map)
                            +EField_Background;
                    #endif

                    OldPosition = Position;
                    double PreviousVelocity = Velocity.getz();
                    UpdateVelocityBoris(SpeciesMass,EField,BField,TimeStep,
                        Velocity,SPEC_CHARGE);

                    #ifdef COLLISIONS
                        double Density = 1.0;
                        double CrossSection = 1.0;
                        double Lambda = 1.0/(Density*CrossSection);
                        //!< Charge exchange collision
                        if( collision_probability(Tau,Velocity.mag3(),Lambda,
                            randnumbers[omp_get_thread_num()]) ){
                                GenerateOrbit(DummyPosition,Velocity,
                                ImpactParameter,ProbUpper,zmin,zmax,0.0,
                                iThermalVel,xThetaPDFMax,
                                randnumbers[omp_get_thread_num()]);
                        }
                    #endif
                    TotalTime+=TimeStep;


                    if( (Velocity.getz()*PreviousVelocity <= 0.0) ){
                        reflections ++;
                    }

                    Position+=TimeStep*Velocity;
                    #ifdef SAVE_APPROACH
                    if( OldPosition.mag3() > Position.mag3()
                        || Position.mag3() < MinPos ){
                        MinPos = Position.mag3();
                    }
                    #endif

                    #ifdef SPHERICAL_INJECTION
                        EdgeCondition = ((Position.getx()*Position.getx()
                            +Position.gety()*Position.gety()
                            +Position.getz()*Position.getz())
                            <= ImpactParameter*ImpactParameter*1.01);
                    #else
                        EdgeCondition = (Position.getz() >= zmin
                            && Position.getz() <= zmax);
                    #endif
                    SphereCondition = (sqrt(Position.getx()*Position.getx()*a1
                        +Position.gety()*Position.gety()*a2
                        +Position.getz()*Position.getz()*a3) > 1.0);

                    #ifdef NO_SPHERE
                        SphereCondition = true;
                    #elif defined DISK
                        SphereCondition
                            = ( (Velocity.getz()*Position.getz() < 0.0)
                            && (sqrt(Position.getx()*Position.getx()*a1
                            +Position.gety()*Position.gety()*a2) < 1.0) );
                    #else
                        SphereCondition
                            = (sqrt(Position.getx()*Position.getx()*a1
                            +Position.gety()*Position.gety()*a2
                            +Position.getz()*Position.getz()*a3) > 1.0);
                    #endif


                    RECORD_TRACK("\n"); RECORD_TRACK(Position);
                    RECORD_TRACK("\t"); RECORD_TRACK(Velocity);
                    RECORD_TRACK("\t");
                    RECORD_TRACK(sqrt(Position.getx()*Position.getx()
                        +Position.gety()*Position.gety()
                        +Position.getz()*Position.getz()));
                }

                CLOSE_TRACK();

                // ************************************************** //

                // ***** PERFORM END OF TRAJECTORY CALCULATIONS      ***** //
                FinalPosition = 0.5*(OldPosition+Position);
                AngularMom = SpeciesMass*(FinalPosition^Velocity);
                #pragma omp critical
                {
                    // In this case it was captured!
                    if( sqrt(Position.getx()*Position.getx()*a1
                        +Position.gety()*Position.gety()*a2
                        +Position.getz()*Position.getz()*a3) < 1.0 ){
                        double AngVelNorm
                            = 5.0*SpeciesMass*MASS/(2.0*DustMass*AngularScalei);
                        AngularVel
                            = (AngVelNorm)*((FinalPosition^Velocity)
                            -(FinalPosition^(TotalAngularVel^FinalPosition)));

                        #ifdef FREE_AXIS
                        unsigned long long n_gen = imax/smax;
                        double TimeScale
                            = n_gen/(2.0*PI*iImpactParameter*iImpactParameter
                            *Radius*Radius*Radius*iDensity*iThermalVel);
                        DustAxis += TotalAngularVel*TimeScale;
                        AxisAVel = TotalAngularVel.getunit();

                        #endif
//                      ADD_F_AMOM(AngularVel*(2.0*DustMass/5.0));
                        TotalInjectedMomCaptured += SpeciesMass*InitialVel;

                        PRINT_FP(fabs(FinalPosition.mag3()-1)); PRINT_FP("\n");
                        TotalAngularVel += AngularVel;
                        TotalAngularMom += AngularMom;

                        j ++; j_ThisSave ++;
                        CapturedCharge += SPEC_CHARGE;

                        PRINT_CHARGE(j)         PRINT_CHARGE("\t")
                        PRINT_CHARGE(PotentialNorm)     PRINT_CHARGE("\t")
                        PRINT_CHARGE(SPEC_CHARGE);  PRINT_CHARGE("\n")
//                      PRINT_AMOM((AngVelNorm)*(FinalPosition^Velocity));
//                      PRINT_AMOM("\t");
//                      PRINT_AMOM((AngVelNorm)
//                      *(FinalPosition^Velocity)*(1.0/Tau)); PRINT_AMOM("\n");
                        ADD_CHARGE()
                        #ifdef VARIABLE_ASCALE
                        TotalAngularVelThisStep += TotalAngularVel;
                        #endif
                        #ifdef VARIABLE_CSCALE
                        TotalChargeInSave += PotentialNorm;
                        #endif
                        //SAVE_LMOM()
                        if(j % num == 0){
                            SAVE_SPOS()
                            SAVE_EPOS()
                            SAVE_AVEL()
                            SAVE_CHA()
                            SAVE_APP()
                            SAVE_AXIS()
                        }
                        // In this case it was trapped!
                    }else if( reflections >= r_max ){
                        TrappedParticles ++;
                        TrappedCharge += SPEC_CHARGE;
                    }else{              // In this case it missed!
                        if(j % num == 0){
                            SAVE_SPOS()
                            SAVE_EPOS()
                            SAVE_APP()
                        }
                        #ifdef SAVE_MOM
                        LinearMomentumSum += SpeciesMass*Velocity;
                        AngularMomentumSum += AngularMom;
                        #endif
                        TotalLostMom += SpeciesMass*Velocity;
                        TotalInjectedMomMissed += SpeciesMass*InitialVel;
                        //SAVE_LMOM()
                        MissedParticles ++;
                        MissedCharge += SPEC_CHARGE;
                    } // END OF if ( Position.mag3() < 1.0 )
                    SAVE_REF()
                    ADD_F_AMOM(SpeciesMass*(Position^Velocity));
                    PRINT_AMOM("FMom = "); PRINT_AMOM(FINAL_AMOM);
                    PRINT_AMOM("\n");
                    C_FINAL_POT(); D_FINAL_POT();
                    PRINT_ENERGY(i); PRINT_ENERGY("\t");
                    PRINT_ENERGY(100*(Velocity.square()
                        /InitialVel.square()-1.0));
                    PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*InitialVel.square()
                        *Radius*Radius/(Tau*Tau));  PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*Velocity.square()
                        *Radius*Radius/(Tau*Tau));  PRINT_ENERGY("\t");
                    PRINT_ENERGY(SPEC_CHARGE*FinalPot);  PRINT_ENERGY("\t");
                    PRINT_ENERGY(SPEC_CHARGE*InitialPot);  PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*InitialVel.square()
                        *Radius*Radius/(Tau*Tau)
                        +SPEC_CHARGE*InitialPot);
                    PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*Velocity.square()
                        *Radius*Radius/(Tau*Tau)+SPEC_CHARGE*FinalPot);
                    PRINT_ENERGY("\t");

                    PRINT_ENERGY((0.5*MASS*SpeciesMass*Velocity.square()
                        *Radius*Radius/(Tau*Tau)+SPEC_CHARGE*FinalPot)/
                        (0.5*MASS*SpeciesMass*InitialVel.square()
                        *Radius*Radius/(Tau*Tau)+SPEC_CHARGE*InitialPot)-1.0);
                    PRINT_ENERGY("\n");
                    TotalNum ++;
                    TotalCharge += SPEC_CHARGE;
                }
            }
        } // END OF PARALLELISED FOR LOOP
        } // END OF PARALLELISED REGION
        // ***** PRINT ANGULAR MOMENTUM AND CHARGE DATA ***** //
        clock_t end = clock();
        double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        RunDataFile << "\n\n***** Save : " << s << " Completed in "
            << elapsd_secs << "s\n\n";

        #ifdef SAVE_CURRENTS
//      Calculate currents for cylindrical geometry with shape factor
        double CyliCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)
            *pow(iImpactParameter,2.0)/(2.0
            *(1-cos(asin(sqrt(1.0/(1+izmax*izmax
            /(iImpactParameter*iImpactParameter))))))
            *(0.5*(TotalNum+TotalCharge))*iDensity);
        double CyleCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)
            *pow(eImpactParameter,2.0)/(2.0
            *(1-cos(asin(sqrt(1.0/(1+ezmax*ezmax
            /(eImpactParameter*eImpactParameter))))))
            *(0.5*(TotalNum-TotalCharge))*eDensity);
//      Calculate currents for cylindrical geometry
        double CylGeoiCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)
            *pow(iImpactParameter,2.0)/((TotalNum+TotalCharge)*iDensity);
        double CylGeoeCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)
            *pow(eImpactParameter,2.0)/((TotalNum-TotalCharge)*eDensity);
//      Calculate currents for Spherical geometry with Bfield
        double SphiCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)
            *pow(iImpactParameter,2.0)/(0.5*(TotalNum+TotalCharge)*iDensity);
        double SpheCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)
            *pow(eImpactParameter,2.0)/(0.5*(TotalNum-TotalCharge)*eDensity);
//      Calculate currents for Spherical geometry without Bfield
        double SphGeoiCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)
            /(0.5*(TotalNum+TotalCharge)*iDensity
            *(1-cos(asin(1.0/iImpactParameter))));
        double SphGeoeCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)
            /(0.5*(TotalNum-TotalCharge)*eDensity
            *(1-cos(asin(1.0/eImpactParameter))));
        #endif

        #ifdef CUSTOM_POTENTIAL
        //Note to self: delete max/min checks
	const std::vector<std::vector<int>> out_of_bounds = this_field_map.get_num_times_outside_domain();
	std::cout<<"Number of instances outside of custom-defined domain:"
		<<"\n\t First Dimension:\n\t\tBelow: "<<out_of_bounds[0][0]
		<< "\n\t\tAbove: "<< out_of_bounds[0][1]
		<<"\n\t Second Dimension:\n\t\tBelow: "<<out_of_bounds[1][0]
		<< "\n\t\tAbove: "<< out_of_bounds[1][1]
		<<"\n\t Third Dimension:\n\t\tBelow: "<<out_of_bounds[2][0]
		<< "\n\t\tAbove: "<< out_of_bounds[2][1]<<std::endl;

        #endif
        #ifdef VARIABLE_CSCALE
        if( j_ThisSave > jmin ){ // Handle no charges captured this save
            MeanChargeDiff = TotalChargeInSave/(j_ThisSave)-MeanChargeSave;
            // If change of sign in Mean Diff, then approaching equilibrium
            if( (MeanChargeDiff < 0 && OldMeanChargeDiff > 0 )
                || (MeanChargeDiff > 0 && OldMeanChargeDiff < 0)  ){
                UPDATE_CSCALE(); // Update charging scale

                // Don't allow charging scale to fall below 1.0e
                if( fabs(ChargeScale) <= 1.0 ){
                    ChargeScale = 1.0;
                    if( jmin == jfin )
                        s = smax;

                    jmin = jfin;
                }
            }
            OldMeanChargeDiff = MeanChargeDiff;
            MeanChargeSave = TotalChargeInSave/(j_ThisSave);
            TotalChargeInSave = 0.0;
            j_ThisSave = 0;
        }
        #endif
        #ifdef VARIABLE_ASCALE
        if( j_ThisSave > jmin ){ // Handle no charges captured this save
            MeanAngularVelDiff
                = TotalAngularVelThisStep*(1.0/j_ThisSave)-MeanAngularVelSave;
            // If change of sign in Mean Diff, then approaching equilibrium
            if( (MeanAngularVelDiff.getz() < 0
                && OldMeanAngularVelDiff.getz() > 0 )
                || (MeanAngularVelDiff.getz() > 0
                && OldMeanAngularVelDiff.getz() < 0 ) ){
                UPDATE_ASCALE();

                // Don't allow angular scale to fall below specified value
                if( fabs(AngularScalee)
                    >= 10.0*MASS*iDensity*sqrt(echarge*iTemp/MASS)
                    *Tau/(Radius) ){
                    AngularScalee
                        = 10.0*MASS*eDensity
                        *sqrt(echarge*iTemp/MASS)*Tau/(Radius);
                    if( jmin == jfin )
                        s = smax;

                    jmin = jfin;
                }
                if( fabs(AngularScalei)
                    >= 10.0*MASS*eDensity*sqrt(echarge*eTemp/MASS)
                    *Tau/(MassRatio*Radius) ){
                    AngularScalei
                        = 10.0*MASS*eDensity
                        *sqrt(echarge*eTemp/MASS)*Tau/(MassRatio*Radius);

                    if( jmin == jfin )
                        s = smax;

                    jmin = jfin;
                }
            }
            OldMeanAngularVelDiff = MeanAngularVelDiff;
            MeanAngularVelSave = TotalAngularVelThisStep*(1.0/j_ThisSave);
            TotalAngularVelThisStep = TotalAngularVelThisStep*0.0;
            j_ThisSave = 0;
            // ************************************************** //
        }
        #endif

        if( j == jmax ){ // We've collected everything, stop saving!
            s = smax;
        }

        SAVE_CHARGE();
        SAVE_MOM();
        SAVE_CURR();
        SAVE_TOT();
        SAVE_LMOM()


        // ************************************************** //

        // If Mean Charge is deviating by less than 0.1%
        // and this is smaller than charge scale length
        //if( fabs(MeanChargeDiff/PotentialNorm) < 0.001
        //  && MeanChargeDiff != 0.0 && ChargeScale/PotentialNorm < 0.01 ){
        //  RunDataFile << "\n\n* Equilibrium Reached! after saves = " << s
        //  << " *";
        //  s = smax;
        //}

        // ***** PRINT CHARGE AND PATH COUNTERS     ***** //
        if( (i_simulated < NumberOfIons || e_simulated < NumberOfElectrons)
            && s == smax ){
            std::cerr << "\nError! Total particle goal was not reached! Data "
                << "may be invalid!";
            RunDataFile << "\n\n* Error! Total particle goal was not reached! "
                << "*";
            RunDataFile << "\n\n*i_sim =  " << i_simulated << "\te_sim = "
                << e_simulated << "* ";
        }

        // ************************************************** //


        // ***** CLOSE DATA FILES           ***** //
        RunDataFile.close();
        CLOSE_CHARGE();
        CLOSE_MOM();
        CLOSE_CURR();
        CLOSE_TOT();
        CLOSE_LMOM()
        CLOSE_AVEL();
        CLOSE_CHA();
        CLOSE_EPOS();
        CLOSE_SPOS();
        CLOSE_APP();
        CLOSE_REF();
        CLOSE_AXIS();

        // ************************************************** //
    }

    // ************************************************** //


    return 0;
}
