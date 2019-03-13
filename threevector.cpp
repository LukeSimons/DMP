/** @file threevector.cpp
 *  @brief Implementation of threevector class for use in physics models
 *
 *  Constructors and printing functions for the threevector class
 *  algebra of vectors with three dimensions.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#include "threevector.h"

threevector::threevector() 
{
    xcoord = 0.0; ycoord = 0.0; zcoord = 0.0;
}

threevector::threevector(double x, double y, double z)
{
    xcoord = x; ycoord = y; zcoord = z;
}

threevector::threevector(double r, double theta, double phiorz, char type)
{
    if( r == 0.0 )
    {
        std::cerr << "Error: threevector Polar parameterised constructor.";
        std::cerr << "Polar coordinate r == 0. theat and phiorz badly defined";
    }
    //! Spherical case
    if(type=='s')
    {
        //! Convert from spherical to cartesian coordinates
        xcoord = r*sin(theta)*cos(phiorz);
        ycoord = r*sin(theta)*sin(phiorz);
        zcoord = r*cos(theta);
    }

    //! Cylindrical case
    if(type=='c')
    {
        //! Convert from cylindrical to cartesian coordinates
        xcoord = r*cos(theta);
        ycoord = r*sin(theta);
        zcoord = phiorz;
    }
    else{
        std::cerr << "Error: Unknown type of polar threevector" << std::endl;
        xcoord = 0.0; ycoord = 0.0; zcoord = 0.0;
    }
}

void threevector::print()
{
    std::cout << "(" << xcoord << ", " << ycoord << ", " << zcoord << ")" 
                << std::endl;
}

void threevector::print(std::ofstream &fout)
{
    fout << "(" << xcoord << ", " << ycoord << ", " << zcoord << ")" << '\t';
}

std::ostream& operator<<(std::ostream& os, const threevector &v)
{
    os << v.getx() << " " << v.gety() << " " << v.getz();
    return os;
}

