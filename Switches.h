/** @file Switches.h
 *  @brief File containing pre-processor directives used to parameterise
 *  the simulation
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs
 */

#define SAVE_TRACKS         //!< Switch: save all particle trajectories
//#define SELF_CONS_CHARGE      //!< Switch  cause incident charges to contribute
//#define VARIABLE_CSCALE       //!< Switch: make particle charges weighted
//#define VARIABLE_ASCALE     //!< Switch: make particle momentum weighted
//#define MAGNETIC_BOTTLE     //!< Switch: Determine collection from potential

#define SAVE_MISSED_MOM     //!< Switch: write missing particles momentum
#define SAVE_ANGULAR_VEL    //!< Switch: write particle charges weighted
#define SAVE_LINEAR_MOM     //!< Switch: write momentum change
#define SAVE_CHARGING         //!< Switch: write the charge collected
#define SAVE_STARTPOS       //!< Switch: write the initial positions
#define SAVE_ENDPOS         //!< Switch: write the final positions
#define SAVE_APPROACH       //!< Switch: write the closest approach pos
#define SAVE_CURRENTS       //!< Switch: write the currents to the sphere
#define SAVE_TOTALS           //!< Switch: write the total currents
#define SAVE_REFLECTS       //!< Switch: write the number of reflections

#define SPHERICAL_INJECTION //!< Switch: inject particles over sphere
//#define POINT_INJECTION     //!< Switch: inject particles at single point
//#define NO_SPHERE           //!< Switch: remove sphere simulation boundary
//#define DISK                //!< Switch: set boundary as flat circle 
//#define COLLISIONS          //!< Switch: Turn neutral Collisions on or off

//#define COULOMB_POTENTIAL   //!< Switch: Coulomb potential for E field
//#define DEBYE_POTENTIAL       //!< Switch: Debye-Huckel potential for E field
//#define PARABOLIC_POTENTIAL //!< Switch: Parabolic potential for E field
#define CUSTOM_POTENTIAL    //!< Switch: Custom potential for E field
//#define BOLTZMANN_DENSITY   //!< Switch: alter injection surface densities

//#define TEST_VELPOSDIST     //!< Print initial velocity and position
//#define TEST_FINALPOS       //!< Print final position of all particles
//#define TEST_CHARGING       //!< Print charge of sphere every particle
//#define TEST_ANGMOM         //!< Print initial and final angular momentum
//#define TEST_COULOMB_ENERGY //!< Print initial and final energies, coulomb
//#define TEST_DEBYE_ENERGY   //!< Print initial and final energies, debye
