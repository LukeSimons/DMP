/** @file Field_Point.cpp
 *  @brief Implementation of field_point class for use in describing Physics fields
 *  field_point Objects contain a value at a point with a three-vector (defined in file threevector.h)
 *
 *  @author Daniel Greenhouse (dg2020@ic.ac.uk (until 01.09.2020), dmgreenhouse@outlook.com (thereafter))
 *      Written immitating the style of main DiMPl author Luke Simons (ls5115@ic.ac.uk)
 *
 *  @bugs No known bugs.
 */

#include "threevector.h" //!< For threevector class
#include "Field_Point.h" //!< For the Field_Point class hereby implemented

/** Constructors */
Field_Point::Field_Point()
{
	value = 0.0; position = threevector();
}

Field_Point::Field_Point( double val, threevector pos)
{
	value = val; position = pos;
}


void Field_Point::add_quadratic_fit(std::vector<std::vector<double>> abc_details){
    // Note to self: requires ability to make not 3D
    _first_d_abc = abc_details[0];
    _second_d_abc = abc_details[1];
    _third_d_abc = abc_details[2];
}

void Field_Point::find_electric_field(){
    // Note to self: requires ability to make 1D, 2D
    //
    // E = -grad (V)

    double E_x = -(2*_first_d_abc[0]*position.getx()+_first_d_abc[1]);
    double E_y = -(2*_second_d_abc[0]*position.gety()+_second_d_abc[1]);
    double E_z = -(2*_third_d_abc[0]*position.getz()+_third_d_abc[1]);
    _electric_field = threevector(E_x, E_y, E_z);
}

threevector Field_Point::find_nearby_E_field(threevector nearby_position){
    // Note to self: requires ability to make 1D, 2D
    //
    // E = -grad (V)
    double E_x = -(2*_first_d_abc[0]*nearby_position.getx()+_first_d_abc[1]);
    double E_y = -(2*_second_d_abc[0]*nearby_position.gety()+_second_d_abc[1]);
    double E_z = -(2*_third_d_abc[0]*nearby_position.getz()+_third_d_abc[1]);
    std::cout<<"First abc:"<<_first_d_abc[0]<<", "<<_first_d_abc[1]<<", "<<_first_d_abc[2]<<std::endl;
    std::cout<<"Second abc:"<<_second_d_abc[0]<<", "<<_second_d_abc[1]<<", "<<_second_d_abc[2]<<std::endl;
    std::cout<<"Third abc:"<<_third_d_abc[0]<<", "<<_third_d_abc[1]<<", "<<_third_d_abc[2]<<std::endl;

    return threevector(E_x, E_y, E_z);
}
