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
	_value = 0.0; _position = threevector();
}

Field_Point::Field_Point( double val, threevector pos)
{
	_value = val; _position = pos;
}

Field_Point::Field_Point( threevector val, threevector pos)
{
	_value = 0.0; // Note to self: change this?
    _electric_field = val;
    _position = pos;
}


void Field_Point::add_quadratic_fit(std::vector<std::vector<double>> abc_details){
    // Note to self: requires ability to make not 3D
    _first_d_abc = abc_details[0];
    _second_d_abc = abc_details[1];
    _third_d_abc = abc_details[2];
}

threevector Field_Point::calculate_electric_field(threevector position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical){
    threevector electric_field;
    // E = -grad (V)
    if (is_cartesian){
        const double E_x = -(2*_first_d_abc[0]*position.getx()+_first_d_abc[1]);
        const double E_y = -(2*_second_d_abc[0]*position.gety()+_second_d_abc[1]);
        const double E_z = -(2*_third_d_abc[0]*position.getz()+_third_d_abc[1]);
        electric_field = threevector(E_x, E_y, E_z);
    } else if (is_spherical){
        const double E_r = -(2*_first_d_abc[0]*position.mag3()+_first_d_abc[1]);
        double E_theta = 0.0; //updated if necessary
        double E_phi = 0.0; //updated if necessary
        if (dimension == 2 || dimension == 3){
            E_theta = -(1/position.mag3())*(2*_second_d_abc[0]*position.gettheta()+_second_d_abc[1]);
            if (dimension == 3) {
                E_phi = -(1/position.mag3()*sin(position.gettheta()))*(2*_third_d_abc[0]*position.getphi()+_third_d_abc[1]);
            }
        }
        electric_field = threevector(E_r, E_theta, E_phi, 's');
    } else if (is_cylindrical){
        const double E_rho = -(2*_first_d_abc[0]*position.getrho()+_first_d_abc[1]);
        double E_z = 0.0; //updated if necessary
        double E_phi = 0.0; //updated if necessary
        if (dimension == 2 || dimension == 3){
            E_z = -(2*_second_d_abc[0]*position.getz()+_second_d_abc[1]);
            if (dimension == 3) {
                E_phi = -(1/position.rho())*(2*_third_d_abc[0]*position.getphi()+_third_d_abc[1]);
            }
        }
        electric_field = threevector(E_r, E_theta, E_phi, 's');
    }
    return electric_field;
}

void Field_Point::find_E_field_at_point(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical){
    _electric_field = calculate_electric_field(_position, dimension, is_cartesian, is_spherical, is_cylindrical);
}

threevector Field_Point::get_nearby_E_field(threevector nearby_position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical){
    // Note to self: this is just a shell for calculate_electric_field

    return calculate_electric_field(nearby_position, dimension, is_cartesian, is_spherical, is_cylindrical);
}
