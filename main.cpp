/** @file main.cpp
 *  @brief DiMPl program to calculate behaviour of conducting spheres in 
 *  plasmas with background flow and magnetic fields
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs
 */


//#define SAVE_TRACKS //!< Switch: save all particle trajectories
//#define SELF_CONS_CHARGE //!< Switch that causes incident charges to contribute
//#define VARIABLE_CSCALE //!< Switch: make particle charges weighted
//#define VARIABLE_ASCALE //!< Switch: make particle momentum weighted

#define SAVE_MISSED_MOM //!< Switch: write to file missing particles momentum
#define SAVE_ANGULAR_VEL //!< Switch: write to file particle charges weighted
#define SAVE_LINEAR_MOM  //!< Switch: write to file momentum change
#define SAVE_CHARGING //!< Switch: write to file the charge collected
#define SAVE_STARTPOS //!< Switch: write to file the initial positions
#define SAVE_ENDPOS //!< Switch: write to file the final positions
//#define SAVE_APPROACH //!< Switch: write to file the closest approach pos
#define SAVE_CURRENTS //!< Switch: write to file the currents to the sphere
#define SAVE_TOTALS //!< Switch: write to file the total currents

#define SPHERICAL_INJECTION //!< Switch: inject particles over sphere
//#define POINT_INJECTION //!< Switch: inject particles at single point
//#define NO_SPHERE //!< Switch: remove inner simulation boundary of sphere

//#define COULOMB_POTENTIAL //!< Switch: use coulomb potential for E field
#define DEBYE_POTENTIAL //!< Switch: use Debye-Huckel potential for E field
//#define BOLTZMANN_DENSITY //!< Switch: alter injection surface densities

//#define TEST_VELPOSDIST //!< Print initial velocity and position distributions
//#define TEST_FINALPOS //!< Print final position of all particles
//#define TEST_CHARGING //!< Print charge of sphere for every particle collected
//#define TEST_ANGMOM //!< Print initial and final angular momentum
//#define TEST_COULOMB_ENERGY //!< Print initial and final energies, coulomb
//#define TEST_DEBYE_ENERGY //!< Print initial and final energies, debye

#include <omp.h>     //!< For parallelisation
#include <iostream>  //!< for std::cout
#include <array>    
#include <random>    //!< for std::normal_distribution<> etc.
#include <fstream>   //!< for std::ofstream
#include <ctime>     //!< for clock()
#include <math.h>    //!< for fabs()
#include <sstream>   //!< for std::stringstream
#include <assert.h>  //!< for assert()
#include <algorithm> //!< for std::min()

#include "Constants.h"  //!< Define Pre-processor Directives and constants
#include "threevector.h"//!< For threevector class
#include "rand_mwts.h"  //!< Functions to generate correct 1-way maxwellian flux

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
    << "\t\tRadius(=1e-6m) DEFAULT,\t\tBy Default, simulate sphere of size 1um in radius\n\n"
    << "\t-s,--spin SPIN\t\t(hz), Specify the initial rotation frequency of Dust grain alligned with magnetic field axis\n"
    << "\t\tSpin(=0.0hz) DEFAULT,\t\tBy Default, simulate sphere with initially zero z angular momentum\n\n"
    << "\t-a1,--semix SEMIX\t\t(arb), Specify the semi-axis for x in dust radii\n"
    << "\t\ta1(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
    << "\t-a2,--semiy SEMIZ\t\t(arb), Specify the semi-axis for y in dust radii\n"
    << "\t\ta2(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
    << "\t-a3,--semiz SEMIY\t\t(arb), Specify the semi-axis for z in dust radii\n"
    << "\t\ta3(=1) DEFAULT,\t\t\tBy Default, simulate perfect sphere\n\n"
    << "\t-d,--density DENSITY\t\t(kgm^-^3), Specify density of Dust grain\n"
    << "\t\tDensity(=19600kgm^-^3) DEFAULT,\tBy Default, Tungsten Density\n\n"
    << "\t-ex,--efieldx EFIELDX\t(double), Specify the x component of a background electric field\n"
    << "\t-ey,--efieldx EFIELDY\t(double), Specify the y component of a background electric field\n"
    << "\t-ez,--efieldx EFIELDZ\t(double), Specify the z component of a background electric field\n"
    << "\t-p,--potential POTENTIAL\t(double), Specify the potential of Dust grain normalised to electron temperature\n"
    << "\t\tPotential(=-2.5eV) DEFAULT,\tBy Default, OML Potential in Ti=Te Hydrogen plasma\n\n"
    << "\t-m,--magfield MAGFIELD\t\t(T), Specify the magnetic field (z direction)\n"
    << "\t\tBMag(=1.0T) DEFAULT,\t\tBy Default, magnetic field is 1.0T upwards in vertical z\n\n"
    << "\t-n,--normalised NORMALISED\t(bool), whether normalisation (following Sonmor & Laframboise) is on or off\n"
    << "\t\tNormalisedVars(=0) DEFAULT,\tFalse, By Default use Tesla and electron volts\n\n"
    << "\t-te,--etemp ETEMP\t\t(eV), Specify the temperature of plasma electrons\n"
    << "\t\teTemp(=1eV) DEFAULT,\n\n"
    << "\t-ne,--edensity EDENSITY\t\t(m^-^3), Specify the plasma electron density\n"
    << "\t\teDensity(=1e18m^-^3) DEFAULT,\tTokamak density\n\n"
    << "\t-ti,--itemp ITEMP\t\t(eV), Specify the temperature of Ions\n"
    << "\t\tiTemp(=1eV) DEFAULT,\n\n"
    << "\t-ni,--idensity IDENSITY\t\t(m^-^3), Specify the plasma ion density\n"
    << "\t\tiDensity(=1e18m^-^3) DEFAULT,\tTokamak density\n\n"
    << "\t-c,--ichance ICHANCE\t\t\t(eV), Specify the probability of generating an Ion\n"
    << "\t\tiChance(=-0.5) DEFAULT,\tFicticious ion generation probability: i.e Self-consistently generate ions & electrons\n\n"
    << "\t-u,--zmaxcoeff ZMAXCOEFF\t(double), The upper limit of simulation domain as number of Coulomb Interaction lengths\n"
    << "\t\tzMaxCoeff(=1.0) DEFAULT,\tNumber of Interaction distances from (0,0,Radius) plane to max of simulation domain\n\n"
    << "\t-l,--zmincoeff ZMINCOEFF\t(double), The lower limit of simulation domain as number of Coulomb Interaction lengths\n"
    << "\t\tzMaxCoeff(=1.0) DEFAULT,\tNumber of Interaction distances from (0,0,Radius) plane to min of simulation domain\n\n"
    << "\t-z,--zboundforce ZBOUNDFORCE\t(double), Force the absolute value of simulation domain upper and lower boundaries\n"
    << "\t\tZBoundForce(=0.0) DEFAULT,\tBy Default, use -u and -l to determine height of injection plane\n\n"
    << "\t-b,--impactpar IMPACTPAR\t(double), Specify the radial limit of simulation domain as number of distances\n"
    << "\t\tImpactPar(=2.0) DEFAULT,\tBy Default, Radial extent of injection is three gyro-radii from centre\n\n"
    << "\t-f,--forceimppar FORCEIMPPAR\t(double), Force the absolute value of simulation radial distance\n"
    << "\t\tForceImpPar(=0.0) DEFAULT,\tBy Default, use -b to determine radial extent of injection plane\n\n"
    << "\t-i,--imax IMAX\t\t\t(int), Specify the number of particles to be launched\n"
    << "\t\timax(=10000) DEFAULT,\t\tBy Default, inject 10,000 particles\n\n"
    << "\t-j,--jmax JMAX\t\t\t(int), Specify the number of particles to be collected (not exceeding imax)\n"
    << "\t\tjmax(=5000) DEFAULT,\t\tBy Default, stop simulation if 5,000 particles are collected\n\n"
    << "\t-rm,--rmax RMAX\t\t\t(int), Specify the number of reflections in z axis before assuming orbit misses\n"
    << "\t\trmax(=15) DEFAULT,\t\tBy Default, stop simulation if Particle reflects in z axis 15 times\n\n"
    << "\t-t,--timestep TIMESTEP\t\t(double), Specify the multiplicative factor for the time step\n"
    << "\t\tTimeStepFactor(=0.0005) DEFAULT,\n\n"
    << "\t-no,--number NUMBER\t\t\t(int), Specify the number of particles to be captured before saving\n"
    << "\t\tNumber(=1000) DEFAULT,\t\tBy Default, store average values after collecting 1000 particles\n\n"
    << "\t-v,--driftvel DRIFTVEL\t\t(m s^-^1), Specify the drift velocity of the plasma\n"
    << "\t\tDriftVel(=0.0) DEFAULT,\t\tBy Default, No drift velocity\n\n"
    << "\t-se,--seed SEED\t\t\t(double), Specify the seed for the random number generator\n"
    << "\t\tSeed(=1.0) DEFAULT,\t\tBy Default, Seed is 1.0\n\n"
    << "\t-sa,--Saves SAVES\t\t\t(int), Specify the number of particles before saving a run\n"
    << "\t\tSaves(=1000) DEFAULT,\t\tBy Default, Save data to Meta datafile once per 1000.0 particles\n\n"
    << "\t-o,--output OUTPUT\t\t(string), Specify the suffix of the output file\n"
    << "\t\tsuffix(='.txt') DEFAULT,\tBy Default, Save data to Data/DiMPl.txt\n\n"
    << std::endl;
    #ifdef VARIABLE_CSCALE
    std::cerr << "\n\nAdditional Options from VARIABLE_CSCALE!\n"
    << "\t-jm,--jmin JMIN\t\t\t(int), Specify the number of particles to be collected before saving charge\n"
    << "\t\tjmin(=20) DEFAULT,\t\tBy Default, re-assess charging scale after collecting 20 charges\n\n"
    << "\t-jf,--jfin JMIN\t\t\t(int), Specify the number of particles to be collected in final save\n"
    << "\t\tjfin(=100) DEFAULT,\t\tBy Default, collect 100 particles in the final save\n\n";
    #endif

}

/** @brief Function to process user command line input
 *  @param argc The number of arguments input by the user
 *  @param argv The character array which contains the users input
 *  @param i A counter indicating the argument which is being processed
 *  @param ss0 The string stream which is used to process command line input
 *  @param Temp The variable of type T which is being filled with user input
 *  @return An integer representing status of
 *
 *  Process user command line input and return status of success of failure
 */
template<typename T> int InputFunction(int &argc, char* argv[], int &i, std::stringstream &ss0, T &Temp){
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

#ifdef SPHERICAL_INJECTION
/** @brief Function to process user command line input
 *  @param v pointer to the three vector of velocity in spherical coordinates
 *  @param phi the angle phi in the spherical coordinate system
 *  @param theta the angle theta in the spherical coordinate system
 *  @param vx pointer to x component of cartesian velocity
 *  @param vy pointer to y component of cartesian velocity
 *  @param vz pointer to z component of cartesian velocity
 *  @author D. M. Thomas (drew.thomas07@imperial.ac.uk)
 *
 *  Convert a velocity in spherical coordinates to a velocity in Cartesian
 *  coordinates, at a position given in spherical coordinates.
 */
void velocity_in_cartesian_coords(double* v, double phi, double theta, double* vx, double* vy, double* vz)
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
    } while (urv1 == 0.0);  //!< `rnd` very occasionally returns zero!
    r = sqrt(-2.0 * (double) logl(urv1));
    theta = 2.0 * PI * urv2;
    *nrv1 = sd * r * cos(theta);
    *nrv2 = sd * r * sin(theta);
}
#endif

/** @brief Generate randomly distributed initial positions and velocities
 *  @param Position reference to three vector with generated position
 *  @param Velocity reference to three vector with generated velocity
 *  @param ImpactParameter radius of injection plane (circle or sphere)
 *  @param ProbUpper probability of injecting particles on upper plane
 *  @param zmin distance to lower injection plane
 *  @param zmax distance to upper injection plane
 *  @param DriftNorm the background flow velocity, shifts distribution
 *  @param ThermalVel the thermal velocity defining width of distribution
 *  @param mt mersenne twister which generates random numbers
 *
 *  Generate the random initial positions and velocities of particles over 
 *  distributions determined by switches. Default is uniformly distributed
 *  positions over circular area with velocities given by rand_mwts(). 
 *  Spherical injection and point injection are also possible
 */
void GenerateOrbit(threevector &Position, threevector &Velocity, const double &ImpactParameter, const double &ProbUpper, const double &zmin, const double zmax, const double DriftNorm, const double ThermalVel, std::mt19937 &mt){

    // ***** DEFINE RANDOM NUMBER GENERATOR ***** //
    //!< See https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    //!< search for "Each component of the velocity vector has a normal distribution"

    std::uniform_real_distribution<double> rad(0.0, 1.0); // IONS
    #ifdef SPHERICAL_INJECTION
        // ***** RANDOMISE POSITION SPHERICALLY ***** //
        double phi_pos   = 2.0*PI*rad(mt);
        double theta_pos = acos(2.0*rad(mt)-1.0);
        Position.setx(ImpactParameter*sin(theta_pos)*cos(phi_pos));
        Position.sety(ImpactParameter*sin(theta_pos)*sin(phi_pos));
        Position.setz(ImpactParameter*cos(theta_pos));
    
    
        // ***** RANDOMISE VELOCITY SPHERICALLY ***** //
/*	POT Style velocity generation
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
        velocity_in_cartesian_coords(v, phi_pos, theta_pos, &(v_temp[0]), &(v_temp[1]), &(v_temp[2]));

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
        Position.setx(0.0);
        Position.sety(0.0);
        Position.setz(zmax);
        
        Velocity.setx(0.0);
        Velocity.sety(0.0);
        Velocity.setz(DriftNorm);
    #else
        // ***** RANDOMISE POSITION CYLINDRICALLY ***** //
        double radial_pos=ImpactParameter*sqrt(rad(mt));
        double theta_pos =2*PI*rad(mt);
        Position.setx(radial_pos*cos(theta_pos));
        Position.sety(radial_pos*sin(theta_pos));
    
        // ***** RANDOMISE VELOCITY CYLINDRICALLY ***** //
        //  Following work from:
        //  Generating equally weighted test particles from the one-way flux of a drifting Maxwellian
        //  T Makkonen, M I Airila & T Kurki-Suonio
        //  http://iopscience.iop.org/article/10.1088/0031-8949/90/1/015204/data

        double invel(0.0);

        if( rad(mt) > ProbUpper ){
            Position.setz(zmax);
            invel = -rand_mwts(-DriftNorm,ThermalVel,mt);   // Generate z-velocity here to determine starting position 
        }else{
            Position.setz(zmin);
            invel = rand_mwts(DriftNorm,ThermalVel,mt); // Generate z-velocity here to determine starting position 
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
static void UpdateVelocityBoris(double MASS, const threevector& Efield, const threevector& BField, double dt, threevector &Velocity, const int SPEC_CHARGE){

    /*magnitude of t, squared*/
    double t_mag2 = (BField*0.5*SPEC_CHARGE*dt).square();

    /*v prime*/
    threevector v_prime = Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt 
                + ((Velocity + Efield*0.5*(SPEC_CHARGE/MASS)*dt)^BField*0.5*SPEC_CHARGE*dt);
    
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
threevector CoulombField(const threevector &Position, double Charge, double COULOMB_NORM){
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
threevector DebyeHuckelField(const threevector &Position, double Charge, double DebyeLength, double COULOMB_NORM){
    if(Charge==0.0) return threevector(0.0,0.0,0.0);
    return COULOMB_NORM*Charge*(1.0/(Position.mag3()))*exp(-(Position.mag3()-1.0)/DebyeLength)*(1.0/Position.mag3()+1.0/DebyeLength)*Position.getunit();
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
    DECLARE_TRACK();
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
    std::ofstream RunDataFile;   //!< Data file for containing the run 
    std::ofstream InputDataFile; //!< Data file for containing the input

    // ************************************************** //


    // ***** DEFINE DUST PARAMETERS         ***** //
    double Radius         = 1e-6;  //!< m, Radius of dust
    double Spin           = 0.0;   //!< hz, Initial rotation rate
    double Density        = 19600; //!< kg m^-^3, Tungsten
    double Potential      = -2.5;  //!< Normalised Potential, 
    double BMag           = 1.0;   //!< Tesla, Magnitude of magnetic field
    double BMagIn         = BMag;  //!< (arb), input magnetic field,in normalised units or Tesla
    bool   NormalisedVars = false; //!< Is magnetic field Normalised according to Sonmor & Laframboise?
    double a1       = 1.0;         //!< Semi-axis for x in dust-grain radii
    double a2       = 1.0;         //!< Semi-axis for y in dust-grain radii 
    double a3       = 1.0;         //!< Semi-axis for z in dust-grain radii

    // ************************************************** //


    // ***** DEFINE PLASMA PARAMETERS       ***** //
    double EFieldx      = 0.0;  //!< V/m, x component of electric field
    double EFieldy      = 0.0;  //!< V/m, x component of electric field
    double EFieldz      = 0.0;  //!< V/m, x component of electric field
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
    #ifdef VARIABLE_CSCALE
    unsigned long long jmin = 20;   //!< Arb, Minimum number of particles collected before summing mean of charge
    unsigned long long jfin = 100;  //!< Arb, Minimum number of particles collected in the final save for averaging
    #endif
    unsigned long long num  = 1000; //!< Arb, Number of particles to be collected before saving
    double TimeStepFactor   = 0.0005; //!< Arb, Multiplicative factor used to determine size of the timestep
    unsigned int Saves(100);    //!<Arb, Number of particles to be collected before saving in a run
    unsigned int reflectionsmax(15); //!< Arb, Number of reflections before rejecting particles


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
        if     ( arg == "--help"    || arg == "-h" ){   show_usage( argv[0]); return 0;         }
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
            InputStatus = InputFunction(argc,argv,i,ss0,NormalisedVars);
        else if( arg == "--etemp"   || arg == "-te")    
            InputStatus = InputFunction(argc,argv,i,ss0,eTemp);
        else if( arg == "--edensity"    || arg == "-ne")    
            InputStatus = InputFunction(argc,argv,i,ss0,eDensity);
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
        #ifdef VARIABLE_CSCALE
        else if( arg == "--jmin"    || arg == "-jm")    
            InputStatus = InputFunction(argc,argv,i,ss0,jmin);
        else if( arg == "--jfin"    || arg == "-jf")    
            InputStatus = InputFunction(argc,argv,i,ss0,jfin);
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
        else if( arg == "--Saves"   || arg == "-sa")    
            InputStatus = InputFunction(argc,argv,i,ss0,Saves);
        else if( arg == "--output"  || arg == "-o" )    
            InputStatus = InputFunction(argc,argv,i,ss0,suffix);
        else{
            sources.push_back(argv[i]);
        }
        if( InputStatus == 1 ){
            std::cerr << "\nError! Input not recognised";
            return 1;
        }
    }
    //!< Input checking
    if( Radius < 0.0 )
        std::cerr << "\nError! Probe Radius is negative\nRadius : " << Radius;
    if( Density < 0.0 )
        std::cerr << "\nError! Probe Density is negative\nDensity : " << Density;
    if( eTemp < 0.0 )
        std::cerr << "\nError! Electron Temperature is negative\neTemp : " << eTemp;
    if( eDensity < 0.0 )
        std::cerr << "\nError! Electron Density is negative\neDensity : " << eDensity;
    if( iTemp < 0.0 )
        std::cerr << "\nError! Ion Temperature is negative\niTemp : " << iTemp;
    if( iDensity < 0.0 )
        std::cerr << "\nError! Ion Density is negative\niDensity : " << iDensity;
    if( iChance < 0.0 )
        std::cout << "\nWarning! Chance of generating Ion is negative, self-consistent flux assumed\niChance : " << iChance;
    if( zMinCoeff < 0.0 )
        std::cerr << "\nError! Lower Vertical Boundary Parameter is negative\nzMinCoeff : " << zMinCoeff;
    if( zMaxCoeff < 0.0 )
        std::cerr << "\nError! Upper Vertical Boundary Parameter is negative\nzMaxCoeff : " << zMaxCoeff;
    if( ZBoundForce < 0.0 )
        std::cerr << "\nError! Force Vertical Boundaries Parameter is negative\nZBoundForce : " << ZBoundForce;
    if( ImpactPar < 0.0 )
        std::cerr << "\nError! Impact Parameter is negative\nImpactPar : " << ImpactPar;
    if( ForceImpPar < 0.0 )
        std::cerr << "\nError! Force Impact Parameter is negative\nForceImpPar : " << ForceImpPar;
    if( imax < jmax )
        std::cerr << "\nError! Total particle goal less than captured particle goal\nimax < jmax : " 
            << imax << " < " << jmax;
    if( jmax < num )
        std::cout << "\nWarning! Save interval less than captured particle goal. No Angular data recorded\nnum < jmax : " 
            << num << " < " << jmax;
    #ifdef VARIABLE_CSCALE
    if( jfin == jmin ){
        std::cout << "\nWarning! Final collected equals minimum collected!\njfin == jmin : " 
            << jfin << " == " << jmin << "\nsetting jfin = jmin+1";
        jfin = jmin+1;
    }
    #endif
    if( Saves > imax ){
        std::cout << "\nwarning! Saves greater than number of simulated particles. Code won't run!\nsetting Saves = imax";
        Saves = imax;
    }

    // ***** WRITE INPUT DATA FILE        ***** //
    //!< Save program configuration to a file with Input
    InputDataFile.open(filename + "_Input" + suffix);
    InputDataFile << "## Input Data File ##\n";
    InputDataFile << "#Input:\nr\ts\ta1\ta2\ta3\td\tp\tm\tn\tte\tne\tti\tni\tc\tu\tl\tz\tb\tf\ti\tj\tt\tno\tv\tse\tsa\to";
    InputDataFile << "\n" << Radius << "\t" << Spin << "\t" << a1 << "\t" << a2 << "\t" << a3 << "\t" << Density  << "\t" << Potential << "\t" << BMagIn << "\t" << NormalisedVars << "\t" << eTemp << "\t" << eDensity<< "\t" << iTemp << "\t" << iDensity << "\t" << iChance << "\t" << zMaxCoeff << "\t" << zMinCoeff << "\t" << ZBoundForce << "\t" << ImpactPar << "\t" << ForceImpPar << "\t" << imax << "\t" << jmax  << "\t" << TimeStepFactor << "\t" << num << "\t" << DriftVel << "\t" << seed << "\t" << Saves << "\t" << suffix;
    InputDataFile.close();
    #ifdef SAVE_TRACKS
    if( imax >= 1000 ){
        std::cout << "\nWarning! Attempting to store more than 1000 tracks!\ni = " << imax;
        std::cout << "\nContinue?...\n";
        std::cin.get();
    }
    #endif

    // ************************************************** //


    //!< If species is positively charged, we assume it's a singly charged ion. Otherwise, singly charged electron
    double MASS = dimplconsts::Mp;  //!< kg, This is the Mass to which quantities are normalised 
    a1 = 1.0/(a1*a1);
    a2 = 1.0/(a2*a2);
    a3 = 1.0/(a3*a3);
    double MassRatio    = sqrt(MASS/Me);
    double DustMass     = (4.0/3.0)*PI*pow(Radius,3)*Density;
    if( NormalisedVars ){   // If we're using S&L normalised units but iChance is undertermined

        if( iChance == 0.0 ){ // If we are simulating only Electrons
                    BMag = sqrt(PI/2.0)*BMagIn*sqrt(Me*eTemp/echarge)/Radius;   // BMag normalised to Electrons
        }else{  // If we are simulating only Ions or otherwise Ions and electrons.
            BMag = sqrt(PI/2.0)*BMagIn*sqrt(MASS*iTemp/echarge)/Radius;   // BMag normalised to Ions
        }
        DriftVel = DriftVel*sqrt(echarge*iTemp/MASS);
    }else{
        BMag = BMagIn;
                Potential = Potential*echarge/(echarge*eTemp);  // Convert from SI Potential to normalised potential
    }

    // ************************************************** //


    // ***** NORMALISATION              ***** //
    // Normalise TIME to the Gyro-Radius of an Ion at B=100T
    // Normalise MASS to Ion Mass
    // Normalise DISTANCE to Dust Radius
    // Normalise CHARGE to fundamental charge
    double MAGNETIC(100);
        double Tau = MASS/(echarge*MAGNETIC);


    double PotentialNorm    = Potential*(eTemp*echarge)*4*PI*epsilon0*Radius/pow(echarge,2.0);  // Normalised Charge,
    double DriftNorm    = DriftVel*Tau/(Radius);
    double DebyeLength  = sqrt((epsilon0*echarge*eTemp)/(eDensity*pow(echarge,2.0)))/Radius;
    double A_Coulomb    = MASS/(4.0*PI*epsilon0*MAGNETIC*MAGNETIC*Radius*Radius*Radius);

    #ifdef VARIABLE_CSCALE
    double ChargeScale  = eTemp*4.0*PI*epsilon0*Radius/(2.0*echarge); // Divide by 2 as we don't want full scale
    #endif
    #ifdef VARIABLE_ASCALE
    double AngularScalei = MASS*iDensity*sqrt(echarge*iTemp/MASS)*Tau*0.01/(Radius);
    DustMass = DustMass*AngularScalei;
    #endif

    // ************************************************** //


    // ***** DEFINE FIELD PARAMETERS        ***** //
    threevector EField_Background(EFieldx,EFieldy,EFieldz);
    threevector Bhat(0.0,0.0,1.0);  // Direction of magnetic field, z dir.

    // ************************************************** //
    

    // ***** DEFINE SIMULATION SPACE        ***** //
    // See : https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    // For reasoning on choice of Velocity vector
    double iThermalVel  = sqrt(echarge*iTemp/MASS)*(Tau/Radius); // Normalised Ion Thermal velocity
    double eThermalVel  = sqrt(echarge*eTemp/Me)*(Tau/Radius);  // Normalised Electron Thermal velocity

    double iRhoTherm = 0.0;                 // Ion gyro-radii are zero by Default
    double eRhoTherm = 0.0;                 // Electron gyro-radii are zero by Default
//  double TimeStepe = TimeStepFactor/fabs(A_Coulomb*MassRatio*MassRatio); // Normalised time step for electrons
//  double TimeStepi = TimeStepFactor/fabs(A_Coulomb);  // Normalised time step for ions
    double TimeStepe = TimeStepFactor;
    double TimeStepi = TimeStepFactor;
    if( BMag == 0.0 && Potential == 0.0 ){
        TimeStepi = 0.01/iThermalVel;
        TimeStepe = 0.01/eThermalVel;
    }
    if( BMag != 0.0 ){                      // If Magnetic field is non-zero
        // Calculate thermal GyroRadius for ions and electrons normalised to dust grain radii
        iRhoTherm   = iThermalVel/(BMag/MAGNETIC); 
        eRhoTherm   = eThermalVel/(pow(MassRatio,2)*BMag/MAGNETIC); 
        
        // Balance electrostatic and magnetic forces to determine impact parameter term

//      iCoulombImpactParameter   = sqrt(fabs(echarge*PotentialNorm
//                      /(4.0*PI*epsilon0*sqrt(echarge*iTemp/Mp)*BMag)))/Radius;
//      eCoulombImpactParameter   = sqrt(fabs(echarge*PotentialNorm
//                      /(4.0*PI*epsilon0*sqrt(echarge*eTemp/Mp)*BMag)))/Radius;

        if( NormalisedVars ){
            iRhoTherm   = 1.0/BMagIn;
            eRhoTherm   = 1.0/BMagIn;
            
            if( iChance != 0.0 ){ // If there is a finite probability of simulating ions
                eRhoTherm       = 1.0/(BMagIn*MassRatio);
            }
        }
        // Time step limited by gyro-motion
        if( TimeStepFactor*eRhoTherm/eThermalVel < TimeStepe ){
            TimeStepe       = TimeStepFactor*eRhoTherm/eThermalVel;   // 1% of Gyro-radius size
        }
        if( TimeStepFactor*iRhoTherm/iThermalVel < TimeStepi ){
                    TimeStepi       = TimeStepFactor*iRhoTherm/iThermalVel;   // 1% of Gyro-radius size
        }
        if( Potential == 0.0 ){
            // Uncharged sphere, time step limited by gyro-motion or thermal velocity
            TimeStepe       = TimeStepFactor*eRhoTherm/eThermalVel;   // 1% of Gyro-radius size
                    TimeStepi       = TimeStepFactor*iRhoTherm/iThermalVel;   // 1% of Gyro-radius size
            if( TimeStepFactor/iThermalVel < TimeStepi ){
                TimeStepi   = TimeStepFactor/iThermalVel;  // Check timestep less than dust radius
            }
            if( TimeStepFactor/eThermalVel < TimeStepe ){
                TimeStepe   = TimeStepFactor/eThermalVel;  // Check timestep less than dust radius
            }
        }
    }

    // Account for non-spherical impact factor, chose larger of two semi axis
    double semiaxisDistortion=1.0;
    if( a1 < a2 ){
        semiaxisDistortion=a1;
    }else{
        semiaxisDistortion=a2;
    }
    double iImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*(iRhoTherm+DebyeLength);
    double eImpactParameter = 1.0/sqrt(semiaxisDistortion)+ImpactPar*(eRhoTherm+DebyeLength);
    if( ForceImpPar > 0.0 ){
        iImpactParameter = 1.0+ForceImpPar;
        eImpactParameter = 1.0+ForceImpPar;
    }

    double iCoulombImpactParameter  = 10.0*fabs(echarge*echarge*PotentialNorm/(2*PI*epsilon0*echarge*iTemp))/Radius; // Balance Coulomb to kinetic energy
    double eCoulombImpactParameter  = 10.0*fabs(echarge*echarge*PotentialNorm/(2*PI*epsilon0*echarge*eTemp))/Radius; // Balance Coulomb to kinetic energy

    double ezmax = 1.05/sqrt(a3)+zMaxCoeff*eCoulombImpactParameter;
    double ezmin = -1.05/sqrt(a3)-zMinCoeff*eCoulombImpactParameter;
    double izmax = 1.0001/sqrt(a3)+zMaxCoeff*iCoulombImpactParameter;
    double izmin = -1.0001/sqrt(a3)-zMinCoeff*iCoulombImpactParameter;
    if( ZBoundForce > 0.0 ){
        ezmax = 1.0/sqrt(a3)+ZBoundForce;
        ezmin = -1.0/sqrt(a3)-ZBoundForce;
        izmax = ezmax;
        izmin = ezmin;
    }

    // ************************************************** //

    // ***** DEFINE PROBABILITY OF ION GENERATION   ***** //
    // Define ratio of flux of electrons to ions
    assert(fabs(ezmax)==fabs(ezmin)); // The min and max heights must match as long as the ProbabilityOfIon is the same for
    assert(fabs(izmax)==fabs(izmin)); // both the top and bottom surfaces.
    #ifdef BOLTZMANN_DENSITY
        #ifdef SPHERICAL_INJECTION
            #ifdef COULOMB_POTENTIAL
            double BoltzmanneDensity = eDensity
                    *exp(PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*eImpactParameter*Radius*echarge*eTemp));
            double BoltzmanniDensity = iDensity
                    *exp(-PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*iImpactParameter*Radius*echarge*iTemp));
            #else
            double BoltzmanneDensity = eDensity
                    *exp(PotentialNorm*echarge*echarge*exp(-(eImpactParameter-1.0)/DebyeLength)
                    /(4.0*PI*epsilon0*eImpactParameter*Radius*echarge*eTemp));
            double BoltzmanniDensity = iDensity
                    *exp(-PotentialNorm*echarge*echarge*exp(-(iImpactParameter-1.0)/DebyeLength)
                    /(4.0*PI*epsilon0*iImpactParameter*Radius*echarge*iTemp));
            #endif
        #else
            #ifdef COULOMB_POTENTIAL
            double BoltzmanneDensity = eDensity
                    *exp(PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
            double BoltzmanniDensity = iDensity
                    *exp(-PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
            #else
            double BoltzmanneDensity = eDensity
                                    *exp(PotentialNorm*echarge*echarge*exp(-(ezmax-1.0)/DebyeLength)
                    /(4.0*PI*epsilon0*ezmax*Radius*echarge*eTemp));
                    double BoltzmanniDensity = iDensity
                                    *exp(-PotentialNorm*echarge*echarge*exp(-(izmax-1.0)/DebyeLength)
                    /(4.0*PI*epsilon0*izmax*Radius*echarge*iTemp));
            #endif
        #endif
    #else
        double BoltzmanneDensity = eDensity;
        double BoltzmanniDensity = iDensity;
    #endif

    double PosFluxi = BoltzmanniDensity*((iThermalVel/sqrt(2.0*PI))*exp(-0.5*DriftNorm*DriftNorm/(iThermalVel*iThermalVel))
            +DriftNorm*0.5*(1.0+erf(DriftNorm/(sqrt(2.0)*iThermalVel))))*pow(iImpactParameter,2);
    double NegFluxi = BoltzmanniDensity*((iThermalVel/sqrt(2.0*PI))*exp(-0.5*DriftNorm*DriftNorm/(iThermalVel*iThermalVel))
            -DriftNorm*0.5*(1.0+erf(-DriftNorm/(sqrt(2.0)*iThermalVel))))*pow(iImpactParameter,2);
    double PosFluxe = BoltzmanneDensity*((eThermalVel/sqrt(2.0*PI))*exp(-0.5*DriftNorm*DriftNorm/(eThermalVel*eThermalVel))
            +DriftNorm*0.5*(1.0+erf(DriftNorm/(sqrt(2.0)*eThermalVel))))*pow(eImpactParameter,2);
    double NegFluxe = BoltzmanneDensity*((eThermalVel/sqrt(2.0*PI))*exp(-0.5*DriftNorm*DriftNorm/(eThermalVel*eThermalVel))
            -DriftNorm*0.5*(1.0+erf(-DriftNorm/(sqrt(2.0)*eThermalVel))))*pow(eImpactParameter,2);

    double TotalFluxProbability = PosFluxi + NegFluxi + PosFluxe + NegFluxe;

    double ProbOfPosFluxe = PosFluxe / TotalFluxProbability;
    double ProbOfNegFluxe = NegFluxe / TotalFluxProbability;
    double ProbOfPosFluxi = PosFluxi / TotalFluxProbability;
    double ProbOfNegFluxi = NegFluxi / TotalFluxProbability;

    double ProbabilityOfIon = ProbOfPosFluxi+ProbOfNegFluxi;
    if( iChance >= 0.0 && iChance <= 1.0 ){
        ProbabilityOfIon = iChance;
        TotalFluxProbability = ProbabilityOfIon*(PosFluxi + NegFluxi) + (1.0-ProbabilityOfIon)*(PosFluxe + NegFluxe);

        ProbOfPosFluxe = (1.0-ProbabilityOfIon)*PosFluxe / TotalFluxProbability;
        ProbOfNegFluxe = (1.0-ProbabilityOfIon)*NegFluxe / TotalFluxProbability;
        ProbOfPosFluxi = ProbabilityOfIon*PosFluxi / TotalFluxProbability;
        ProbOfNegFluxi = ProbabilityOfIon*NegFluxi / TotalFluxProbability;
    }

    // ************************************************** //


    // ***** SEED RANDOM NUMBER GENERATOR IN THREADS***** //
    std::random_device rd;      // Create Random Device
    std::vector<std::mt19937> randnumbers;
    std::uniform_real_distribution<double> rad(0, 1); // IONS
    for(int p = 0; p < omp_get_max_threads(); p ++){
        randnumbers.push_back(std::mt19937(seed+p));
    }

    // ************************************************** //


    // ***** OPEN DATA FILE WITH HEADER         ***** //
    time_t now = time(0);       // Get the time of simulation
    char * dt = ctime(&now);

    OPEN_CHARGE();  HEAD_CHARGE();
    OPEN_MOM(); HEAD_MOM();
    OPEN_AVEL();    HEAD_AVEL();
    OPEN_LMOM();    HEAD_LMOM();
    OPEN_CHA(); HEAD_CHA();
    OPEN_SPOS();    HEAD_SPOS();
    OPEN_EPOS();    HEAD_EPOS();
    OPEN_APP(); HEAD_APP()
    OPEN_CURR();    HEAD_CURR();
    OPEN_TOT(); HEAD_TOT();

    unsigned long long NumberOfIons = ProbabilityOfIon*imax;
    unsigned long long NumberOfElectrons = imax-NumberOfIons;
    
    RunDataFile.open(filename + suffix);
    RunDataFile << "## Run Data File ##\n";
    RunDataFile << "#Date: " << dt;
    RunDataFile << "#Input:\t\tValue\n\nimax (arb #):\t\t"<<imax<<"\njmax (arb #):\t\t"<<jmax<<"\nIon # (arb #):\t\t" << NumberOfIons<<"\nElec # (arb #):\t\t"<<NumberOfElectrons<<"\nProbOfIon (arb):\t"<<ProbabilityOfIon<<"\n\nElectron Gyro (1/Radius):\t"<<eRhoTherm<<"\nElectron Temp (eV):\t\t"<<eTemp<<"\nElec Density (m^-^3):\t\t"<<BoltzmanneDensity<<"\nElectron IP (1/Radius):\t\t"<<eImpactParameter<<"\nElectron zmax (1/Radius):\t"<<ezmax<<"\nElectron zmin (1/Radius):\t"<<ezmin<<"\nElecs Timestep (Tau):\t\t"<<TimeStepe<<"\nProbOfPosFluxe (arb):\t\t"<<ProbOfPosFluxe<<"\nProbOfNegFluxe (arb):\t\t"<<ProbOfNegFluxe<<"\n\nIon Gyro (1/Radius):\t"<<iRhoTherm<<"\nIon Temp (eV):\t\t"<<iTemp<<"\nIon Density (m^-^3):\t"<<BoltzmanniDensity<<"\nIon IP (1/Radius):\t"<<iImpactParameter<<"\nIon zmax (1/Radius):\t"<<izmax<<"\nIon zmin (1/Radius):\t"<<izmin<<"\nIon Timestep (Tau):\t"<<TimeStepi<<"\nProbOfPosFluxi (arb):\t"<<ProbOfPosFluxi<<"\nProbOfNegFluxi (arb):\t"<<ProbOfNegFluxi<<"\n\nRadius (m):\t\t"<<Radius<<"\nSpin (1/Tau):\t\t"<<Spin*Tau<<"\na1 (1/Radius):\t\t"<<(1.0/sqrt(a1))<<"\na2 (1/Radius):\t\t"<<(1.0/sqrt(a2))<<"\na3 (1/Radius):\t\t"<<(1.0/sqrt(a3))<<"\nDensity (kg m^-^3):\t"<<Density<<"\nCharge (1/echarge):\t\t"<<PotentialNorm<<"\nB Field (T or Radius/GyroRad):\t"<<BMag<<"\nDebyeLength (1/Radius):\t\t"<<DebyeLength <<"\nDrift Norm (Radius/Tau):\t"<<DriftNorm<<"\nTime Norm [Tau] (s):\t\t"<<Tau<<"\n\n"<<"RNG Seed (arb):\t\t"<<seed<<"\nOMP_THREADS (arb):\t"<<omp_get_max_threads()<<"\n\n";
    #ifdef SPHERICAL_INJECTION
        RunDataFile << "* SPHERICAL INJECTION *\n";
    #else 
        RunDataFile << "* CYLINDRICAL INJECTION *\n";
    #endif
    #ifdef NO_SPHERE
        RunDataFile << "* NO SPHERE *\n";
    #endif
    #ifdef DEBYE_POTENTIAL
        RunDataFile << "* DEBYE HUCKEL POTENTIAL *\n\n";
    #endif
    #ifdef COULOMB_POTENTIAL
        RunDataFile << "* COULOMB POTENTIAL *\n\n";
    #endif
    #ifdef BOLTZMANN_DENSITY
        RunDataFile << "* BOLTZMANN DENSITY *\n\n";
    #endif  
    // ************************************************** //


    // ***** BEGIN LOOP OVER PARTICLE ORBITS    ***** //
    threevector TotalAngularVel(0.0,0.0,Spin*Tau);

    threevector MeanAngularVel(0.0,0.0,Spin*Tau);
    threevector TotalAngularVelThisStep(0.0,0.0,Spin*Tau);
    threevector MeanAngularVelDiff(0.0,0.0,0.0);

    threevector TotalAngularMom(0.0,0.0,0.0);
    threevector TotalInjectedMomCaptured(0.0,0.0,0.0);
	threevector TotalInjectedMomMissed(0.0,0.0,0.0);
    threevector TotalLostMom(0.0,0.0,0.0);

    #ifdef VARIABLE_CSCALE
    double MeanChargeSave = PotentialNorm;  // Initialise the mean charge as starting charge
    double TotalChargeInSave=PotentialNorm;// Initialise the charge collected in this save
    double MeanChargeDiff(0.0);         // Initial Mean charge diff is zero
    double OldMeanChargeDiff(0.0);      // Initial Mean charge diff is zero
    #endif  

    #ifdef SAVE_MOM
    threevector LinearMomentumSum(0.0,0.0,0.0);
    threevector AngularMomentumSum(0.0,0.0,0.0);
    #endif
    
    #ifdef TEST_ANGMOM
    threevector INITIAL_AMOM(0.0,0.0,0.0);
    threevector FINAL_AMOM(0.0,0.0,0.0);
    #endif

    unsigned long long j(0), i(0), RegeneratedParticles(0), TrappedParticles(0), MissedParticles(0), TotalNum(0);
    unsigned long long j_ThisSave(0), e_simulated(0), i_simulated(0);
    long long CapturedCharge(0), RegeneratedCharge(0), TrappedCharge(0), MissedCharge(0), TotalCharge(0);

    unsigned long long smax = imax / Saves;

    for( unsigned int s=1; s <= smax; s ++){
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

        // ************************************************** //

        #pragma omp parallel shared(TotalAngularVel,TotalAngularMom,TotalInjectedMomCaptured,TotalInjectedMomMissed,TotalLostMom,j) PRIVATE_FILES()
        {
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
        unsigned int reflections(0);

        double BMagNorm;
        double ImpactParameter;
        double ThermalVel;
        double zmax;   
        double zmin;   
        double TimeStep;
        double SpeciesMass;
        double ProbUpper;

        int SPEC_CHARGE=-1;

        bool EdgeCondition = true;
        bool SphereCondition = true;

        #pragma omp for
        for( i=(IMAX-imax/smax); i < IMAX; i ++){   // Loop over maximum number of particles to generate
            if( j <= jmax ){    // Loop until we reach a certain number of particles jmax
    
                //std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();
    
                // ***** DETERMINE IF IT'S AN ELECTRON OR ION ***** //
                // MassRatio is squared because of absence of mass in Boris solver
                #pragma omp critical
                {
                if( rad(randnumbers[omp_get_thread_num()]) < ProbabilityOfIon && i_simulated < NumberOfIons ){ 
                    // If this is the case, we need to generate an ion
                    BMagNorm = BMag/MAGNETIC;
                    ImpactParameter=iImpactParameter;
                    ThermalVel=iThermalVel;
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
                        zmax= ezmax;        // Top of Simulation Domain, in Dust Radii
                        zmin= ezmin;        // Top of Simulation Domain, in Dust Radii
                        TimeStep = TimeStepe;
                        SpeciesMass = 1.0/pow(MassRatio,2);
                        ProbUpper = ProbOfPosFluxe/(ProbOfPosFluxe+ProbOfNegFluxe);
                        SPEC_CHARGE=-1;
                        e_simulated = e_simulated + 1;
                    }else if( i_simulated < NumberOfIons ){
                        BMagNorm = BMag/MAGNETIC;
                        ImpactParameter=iImpactParameter;
                        ThermalVel=iThermalVel;
                        zmax    = izmax; 
                        zmin    = izmin ;
                        TimeStep = TimeStepi;
                        SpeciesMass = 1.0;
                        ProbUpper = ProbOfPosFluxi/(ProbOfPosFluxi+ProbOfNegFluxi);
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
                // ***** GENERATE AN ORBIT ***** //

                }
                GenerateOrbit(Position,Velocity,ImpactParameter,ProbUpper,zmin,zmax,DriftNorm,ThermalVel,randnumbers[omp_get_thread_num()]);
    
                InitialPos = Position;
                InitialVel = Velocity;
                // ************************************************** //
    
                // ***** TESTING AREA               ***** //
                // ***** ANGULAR-MOMENTUM TEST          ***** //
                #ifdef TEST_ANGMOM
                #pragma omp critical 
                {
                    ADD_I_AMOM(SpeciesMass*(Position^Velocity));        // For Angular Momentum Calculations
                    PRINT_AMOM("IMom = "); PRINT_AMOM(INITIAL_AMOM); PRINT_AMOM("\n");
                }
                #endif
                // ***** VELOCITY-POSITION DISTRIBUTION TEST    ***** //
                #ifdef TEST_VELPOSDIST
                #pragma omp critical
                {
                    PRINT_VPD(Position); PRINT_VPD("\t");   // For debugging look at initial positions
                    PRINT_VPD(Velocity); PRINT_VPD("\t");   // For debugging look at initial velocities
                    PRINT_VPD(Velocity*(Position.getunit())); PRINT_VPD("\t");  
                    PRINT_VPD( sqrt(pow(Velocity.getx(),2)+pow(Velocity.gety(),2))*SpeciesMass); 
                    PRINT_VPD("\n");    // For debugging look at gyro-radii
                }
                #endif
                // ***** ENERGY TEST: MEASURE INITIAL ENERGY    ***** //
                C_INITIAL_VEL(); D_INITIAL_VEL();   // For energy calculations
                C_INITIAL_POT(); D_INITIAL_POT();   // For energy calculations
                #ifdef SAVE_APPROACH
                    double MinPos = sqrt(zmax*zmax+ImpactParameter*ImpactParameter);
                #endif


                // ***** RECORD TRACK DATA, DEBUG AND TEST  ***** //
                OPEN_TRACK(filename + "_Track_" + std::to_string(i) + ".txt");
                RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
    
                // ************************************************** //
    
    
                // ***** TAKE INITIAL HALF STEP BACKWARDS ***** //
                // Calculate Electric Field
                #ifdef DEBYE_POTENTIAL
                    EField = DebyeHuckelField(Position,PotentialNorm,DebyeLength,A_Coulomb)+EField_Background;
                #endif
                #ifdef COULOMB_POTENTIAL
                    EField = CoulombField(Position,PotentialNorm,A_Coulomb) +EField_Background;
                #endif
                UpdateVelocityBoris(SpeciesMass,EField,BField,-0.5*TimeStep,Velocity,SPEC_CHARGE);  
    
                // ************************************************** //
    
    
                // ***** DO PARTICLE PATH INTEGRATION       ***** //
                OldPosition.setx(0.0);
                OldPosition.sety(0.0);
                OldPosition.setz(0.0);

                #ifdef SPHERICAL_INJECTION
                    EdgeCondition = ((Position.getx()*Position.getx()+Position.gety()*Position.gety()+Position.getz()*Position.getz()) <= ImpactParameter*ImpactParameter*1.01);
                #else
                    EdgeCondition = (Position.getz() >= zmin && Position.getz() <= zmax);
                #endif
                
                #ifdef NO_SPHERE
                    SphereCondition = true;
                #else
                    SphereCondition = (sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2+Position.getz()*Position.getz()*a3) > 1.0);
                #endif

                // While we don't exceed a specified number of reflections to catch trapped orbits AND  
                // while the particle is not inside the sphere and not outside the simulation domain
                reflections=0;
                unsigned int r_max = reflectionsmax;
                if( PotentialNorm*SPEC_CHARGE > 0.0 ){ // In this case, potential is repulsive
                    r_max = 1;
/*                  // Compare vertical kinetic energy to potential. If it's smaller, reject orbit
                    // immediately before simulating
                    #ifdef DEBYE_POTENTIAL
                    if( fabs((PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*Radius))
                                               *(1.0-exp((zmax/DebyeLength)-(1.0/DebyeLength))/zmax))
                        > 0.5*SpeciesMass*Mp*Velocity.getz()*Velocity.getz()*Radius*Radius/(Tau*Tau) ){
                            EdgeCondition = false;
                        }
                    #endif
                    #ifdef COULOMB_POTENTIAL
                        if( fabs((PotentialNorm*echarge*echarge/(4.0*PI*epsilon0*Radius))*(1.0-(1.0/zmax)))
                        > (0.5*SpeciesMass*Mp*Velocity.getz()*Velocity.getz()*Radius*Radius/(Tau*Tau)) ){
                            EdgeCondition = false;
                        }
                    #endif*/
                }

                while( SphereCondition && EdgeCondition && reflections < r_max ){
                    #ifdef DEBYE_POTENTIAL
                        EField = DebyeHuckelField(Position,PotentialNorm,DebyeLength,A_Coulomb)+EField_Background;
                    #endif
                    #ifdef COULOMB_POTENTIAL
                        EField = CoulombField(Position,PotentialNorm,A_Coulomb)+EField_Background;
                    #endif
                    OldPosition = Position;
                    double PreviousVelocity = Velocity.getz();
                    UpdateVelocityBoris(SpeciesMass,EField,BField,TimeStep,Velocity,SPEC_CHARGE);
                    

                    if( (Velocity.getz()*PreviousVelocity <= 0.0) ){
                        reflections ++;
                    }

                    Position+=TimeStep*Velocity;
                    #ifdef SAVE_APPROACH
                    if( OldPosition.mag3() > Position.mag3() || Position.mag3() < MinPos ){
                        MinPos = Position.mag3();
                    }
                    #endif

                    #ifdef SPHERICAL_INJECTION
                        EdgeCondition = ((Position.getx()*Position.getx()+Position.gety()*Position.gety()+Position.getz()*Position.getz()) <= ImpactParameter*ImpactParameter*1.01);
                    #else 
                        EdgeCondition = (Position.getz() >= zmin && Position.getz() <= zmax);
                    #endif
                    SphereCondition = (sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2+Position.getz()*Position.getz()*a3) > 1.0);
                
                    #ifdef NO_SPHERE
                        SphereCondition = true;
                    #else
                        SphereCondition = (sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2+Position.getz()*Position.getz()*a3) > 1.0);
                    #endif


                    RECORD_TRACK("\n");RECORD_TRACK(Position);RECORD_TRACK("\t");RECORD_TRACK(Velocity);
                    RECORD_TRACK("\t");
                    RECORD_TRACK(sqrt(Position.getx()*Position.getx()
                            +Position.gety()*Position.gety()+Position.getz()*Position.getz()));
                }   
                CLOSE_TRACK();
                // ************************************************** //
    
    
                // ***** PERFORM MOMENTUM CALCULATIONS      ***** //
                FinalPosition = 0.5*(OldPosition+Position);
                AngularMom = SpeciesMass*(FinalPosition^Velocity);      
                #pragma omp critical
                {
                    if( sqrt(Position.getx()*Position.getx()*a1+Position.gety()*Position.gety()*a2+Position.getz()*Position.getz()*a3) < 1.0 ){ // In this case it was captured!
                        double AngVelNorm = 5.0*SpeciesMass*MASS/(2.0*DustMass);
                        AngularVel = (AngVelNorm)*
                        ((FinalPosition^Velocity)-(FinalPosition^(TotalAngularVel^FinalPosition)));
//                      ADD_F_AMOM(AngularVel*(2.0*DustMass/5.0));
						TotalInjectedMomCaptured += SpeciesMass*InitialVel;

                        PRINT_FP(fabs(FinalPosition.mag3()-1)); PRINT_FP("\n");
                        TotalAngularVel += AngularVel;
                        TotalAngularVelThisStep += TotalAngularVel;
                        
                        TotalAngularMom += AngularMom;

                        j ++; j_ThisSave ++;
                        CapturedCharge += SPEC_CHARGE;

                        PRINT_CHARGE(j)         PRINT_CHARGE("\t")
                        PRINT_CHARGE(PotentialNorm)     PRINT_CHARGE("\t")
                        PRINT_CHARGE(SPEC_CHARGE);  PRINT_CHARGE("\n")
//                      PRINT_AMOM((AngVelNorm)*(FinalPosition^Velocity)); PRINT_AMOM("\t");
//                      PRINT_AMOM((AngVelNorm)*(FinalPosition^Velocity)*(1.0/Tau)); PRINT_AMOM("\n");
                        ADD_CHARGE()
                        #ifdef VARIABLE_CSCALE
                        TotalChargeInSave += PotentialNorm;
                        #endif
                        SAVE_SPOS()
                        SAVE_EPOS()
                        //SAVE_LMOM()
                        if(j % num == 0){
                            SAVE_AVEL()
                            SAVE_CHA()
                            SAVE_APP()
                        }
                    }else if( reflections >= r_max ){        // In this case it was trapped!
                        TrappedParticles ++; 
                        TrappedCharge += SPEC_CHARGE;
                                        }else{              // In this case it missed!
                        SAVE_SPOS()
                        SAVE_EPOS()
                        SAVE_APP()
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

                    ADD_F_AMOM(SpeciesMass*(Position^Velocity));
                    PRINT_AMOM("FMom = "); PRINT_AMOM(FINAL_AMOM); PRINT_AMOM("\n");
                    C_FINAL_POT(); D_FINAL_POT();
                    PRINT_ENERGY(i); PRINT_ENERGY("\t"); 
                    PRINT_ENERGY(100*(Velocity.square()/InitialVel.square()-1.0));  PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*InitialVel.square()*Radius*Radius/(Tau*Tau));  PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*Velocity.square()*Radius*Radius/(Tau*Tau));  PRINT_ENERGY("\t");
                    PRINT_ENERGY(SPEC_CHARGE*FinalPot);  PRINT_ENERGY("\t");
                    PRINT_ENERGY(SPEC_CHARGE*InitialPot);  PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*InitialVel.square()*Radius*Radius/(Tau*Tau)
                            +SPEC_CHARGE*InitialPot);  PRINT_ENERGY("\t");
                    PRINT_ENERGY(0.5*MASS*SpeciesMass*Velocity.square()*Radius*Radius/(Tau*Tau)
                            +SPEC_CHARGE*FinalPot);  PRINT_ENERGY("\t");

                    PRINT_ENERGY((0.5*MASS*SpeciesMass*Velocity.square()*Radius*Radius/(Tau*Tau)
                            +SPEC_CHARGE*FinalPot)/
                            (0.5*MASS*SpeciesMass*InitialVel.square()*Radius*Radius/(Tau*Tau)
                            +SPEC_CHARGE*InitialPot) - 1.0);  
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
        RunDataFile << "\n\n***** Save : " << s << " Completed in " << elapsd_secs << "s\n\n";
    
//      Calculate currents for cylindrical geometry with shape factor
        double CyliCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
        /(2.0*(1-cos(asin(sqrt(1.0/(1+izmax*izmax/(iImpactParameter*iImpactParameter))))))*(0.5*(TotalNum+TotalCharge))*iDensity);
        double CyleCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
        /(2.0*(1-cos(asin(sqrt(1.0/(1+ezmax*ezmax/(eImpactParameter*eImpactParameter))))))*(0.5*(TotalNum-TotalCharge))*eDensity);
//      Calculate currents for cylindrical geometry
        double CylGeoiCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
                    /((TotalNum+TotalCharge)*iDensity);
        double CylGeoeCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
                /((TotalNum-TotalCharge)*eDensity);
//      Calculate currents for Spherical geometry with Bfield
        double SphiCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)*pow(iImpactParameter,2.0)
                    /(0.5*(TotalNum+TotalCharge)*iDensity);
        double SpheCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)*pow(eImpactParameter,2.0)
                /(0.5*(TotalNum-TotalCharge)*eDensity);
//      Calculate currents for Spherical geometry without Bfield
        double SphGeoiCurr = 0.5*BoltzmanniDensity*(j+CapturedCharge)
                    /(0.5*(TotalNum+TotalCharge)*iDensity*(1-cos(asin(1.0/iImpactParameter))));
        double SphGeoeCurr = 0.5*BoltzmanneDensity*(j-CapturedCharge)
                /(0.5*(TotalNum-TotalCharge)*eDensity*(1-cos(asin(1.0/eImpactParameter))));

        #ifdef VARIABLE_CSCALE
        if( j_ThisSave > jmin ){ // Handle no charges captured this save
            MeanChargeDiff = TotalChargeInSave/(j_ThisSave)-MeanChargeSave;
            // If change of sign in Mean Diff, then approaching equilibrium
            if( (MeanChargeDiff < 0 && OldMeanChargeDiff > 0 )
                || (MeanChargeDiff > 0 && OldMeanChargeDiff < 0)  ){
                UPDATE_CSCALE(); // Update charging scale
                
                // Don't allow charging scale to fall below 1.0e
                if( fabs(ChargeScale) <= 2.0 ){
                    ChargeScale = 2.0;
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
        #ifdef VARIABlE_ASCALE
        if( j_ThisSave > jmin ){ // Handle no charges captured this save
            MeanAngularVelDiff = TotalAngularVelThisStep*(1.0/j_ThisSave)-MeanAngularVelSave;

            if( (MeanAngularVelDiff.getz() < 0 && 
                TotalAngularVelThisStep.getz()*(1.0/j_ThisSave)-MeanAngularVelSave.getz() > 0 )
                || (MeanAngularVelDiff.getz() > 0 && 
                TotalAngularVelThisStep.getz()*(1.0/j_ThisSave)-MeanAngularVelSave.getz() < 0 ) ){
                UPDATE_ASCALE();
            }

            MeanAngularVelSave = TotalAngularVelThisStep*(1.0/j_ThisSave);
            TotalAngularVelThisStep = TotalAngularVelThisStep*0.0;
            j_ThisSave = 0;
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
        //if( fabs(MeanChargeDiff/PotentialNorm) < 0.001 && MeanChargeDiff != 0.0 
        //  && ChargeScale/PotentialNorm < 0.01 ){ 
        //  RunDataFile << "\n\n* Equilibrium Reached! after saves = " << s << " *";
        //  s = smax;
        //}
    
        // ***** PRINT CHARGE AND PATH COUNTERS     ***** //
        if( (i_simulated < NumberOfIons || e_simulated < NumberOfElectrons) && s == smax ){
            std::cerr << "\nError! Total particle goal was not reached! Data may be invalid!";
            RunDataFile << "\n\n* Error! Total particle goal was not reached! *";
            RunDataFile << "\n\n*i_sim =  " << i_simulated << "\te_sim = " << e_simulated << "* ";
        }
    
        // ************************************************** //
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
        // ***** CLOSE DATA FILES           ***** //
    }

    // ************************************************** //


    return 0;
}
