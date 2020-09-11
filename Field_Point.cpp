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


Field_Point::Field_Point()
{
    const threevector default_threevec = threevector(); // Initialise a default threevector
	_value = 0.0;
    _position = default_threevec;
}

Field_Point::Field_Point( const double val, threevector pos)
{
	_value = val;
    _position = pos;
}

Field_Point::Field_Point( double one_d_val, double two_d_val, double three_d_val, threevector pos, bool is_cartesian, bool is_spherical, bool is_cylindrical)
{
	_value = 0.0;
    _one_d_val = one_d_val;
    _two_d_val = two_d_val;
    _three_d_val = three_d_val;
    _position = pos;
}


void Field_Point::add_quadratic_fit(std::vector<std::vector<double>> abc_details){
    _first_d_abc = abc_details[0];
    _second_d_abc = abc_details[1];
    _third_d_abc = abc_details[2];
}

threevector Field_Point::calculate_electric_field(threevector position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical){
    double one_d_val; // first component of the electric field
    double two_d_val; // second component of the electric field
    double three_d_val; // third component of the electric field
    // E = -grad (V)
    if (is_cartesian){
        // According to x direction: V = ax^2 + bx + c
        one_d_val = -(2*_first_d_abc[0]*position.getx()+_first_d_abc[1]);
        // According to y direction: V = ay^2 + by + c
        two_d_val = -(2*_second_d_abc[0]*position.gety()+_second_d_abc[1]);
        // According to z direction: V = az^2 + bz + c
        three_d_val = -(2*_third_d_abc[0]*position.getz()+_third_d_abc[1]);
    } else if (is_spherical){
        // Find position at which the electric field is to be calculated in spherical coordinates
        const double pos_r = position.mag3();
        const double pos_theta = position.gettheta();
        const double pos_phi = position.getphi();

        // Declare electric field in spherical coordinates (these later become one_d_val, two_d_val, three_d_val)
        // According to radial direction: V = ar^2 + br + c
        const double E_r = -(2*_first_d_abc[0]*pos_r+_first_d_abc[1]);
        double E_theta;
        double E_phi;
        if (dimension==1){
            // According to theta direction, no variation.
            E_theta = 0.0;
            // According to phi direction, no variation.
            E_phi = 0.0;
        } else if (dimension == 2){
            // According to theta direction: V = a theta^2 + b theta + c
            // Necessary Jacobian accounted for.
            E_theta = -(1/pos_r)*(2*_second_d_abc[0]*pos_theta+_second_d_abc[1]);
            // According to phi direction, no variation.
            E_phi = 0.0;
        } else if (dimension == 3){
            // According to theta direction: V = a theta^2 + b theta + c
            // Necessary Jacobian accounted for.
            E_theta = -(1/pos_r)*(2*_second_d_abc[0]*pos_theta+_second_d_abc[1]);

            double modified_pos_theta = pos_theta; //if theta is 0 or pi, shift slightly to avoid infinity
            const double delta = 1e-15;
            if (pos_theta==0.0){
                modified_pos_theta+=delta;
            } else if (pos_theta>=PI){
                modified_pos_theta-=delta;
            }
            // According to phi direction: V = a phi^2 + b phi + c
            // Necessary Jacobian accounted for.
            E_phi = -(1/(sin(modified_pos_theta)*pos_r))*(2*_third_d_abc[0]*pos_phi+_third_d_abc[1]);
        }
        one_d_val = E_r;
        two_d_val = E_theta;
        three_d_val = E_phi;
    } else if (is_cylindrical){
        // Find position at which the electric field is to be calculated in Cylindrical coordinates
        const double pos_rho = position.getrho();
        const double pos_z = position.getz();
        const double pos_phi = position.getphi();

        // Declare the gradient of the potential in cylindrical coordinates along each dimension (the negatives of these later become the electric field component values one_d_val, two_d_val, three_d_val)
        // According to rho direction: V = a rho^2 + b rho + c
        double grad_rho = (2*_first_d_abc[0]*pos_rho+_first_d_abc[1]);
        double grad_phi;
        double grad_z;

        if (dimension==1){
            // According to z direction, no variation.
            grad_z = 0.0;
            // According to phi direction, no variation.
            grad_phi = 0.0;
        } else if (dimension == 2){
            // According to z direction: V = az^2 + bz + c
            grad_z = (2*_second_d_abc[0]*pos_z+_second_d_abc[1]);
            // According to phi direction, no variation.
            grad_phi = 0.0;
        } else if (dimension == 3){
            // According to z direction: V = az^2 + bz + c
            grad_z = (2*_second_d_abc[0]*pos_z+_second_d_abc[1]);
            // According to phi direction: V = a phi^2 + b phi + c
            // Necessary Jacobian accounted for.
            grad_phi = (1/pos_rho)*(2*_third_d_abc[0]*pos_phi+_third_d_abc[1]);
        }
        // Convert gradients to the electric field components
        one_d_val = -grad_rho;
        two_d_val = -grad_z;
        three_d_val = -grad_phi;
    }
    // Ensure returned electric field coorectly fits DiMPl's Cartesian understanding
    return convert_electric_field(one_d_val, two_d_val, three_d_val, position, is_cartesian, is_spherical, is_cylindrical);
}

threevector Field_Point::convert_electric_field(double one_d_val, double two_d_val, double three_d_val, threevector position, bool is_cartesian, bool is_spherical, bool is_cylindrical){
    // Declare Cartesian electric field components.
    double E_x;
    double E_y;
    double E_z;
    if (is_cartesian){
        // Electric Field components already in Cartesian form
        E_x = one_d_val;
        E_y = two_d_val;
        E_z = three_d_val;
    } else if (is_spherical){
        // Find positions in spherical coordinates
        const double pos_theta = position.gettheta();
        const double pos_phi = position.getphi();

        // Declare electric field in spherical coordinates
        const double E_r = one_d_val;
        const double E_theta = two_d_val;
        const double E_phi = three_d_val;

        // For efficiency, find these once and store
        const double sin_pos_theta = sin(pos_theta);
        const double cos_pos_theta = cos(pos_theta);
        const double sin_pos_phi = sin(pos_phi);
        const double cos_pos_phi = cos(pos_phi);

        // Transform to cartesian basis
        E_x = (sin_pos_theta*cos_pos_phi*E_r)+(cos_pos_theta*cos_pos_phi*E_theta)+(-sin_pos_phi*E_phi);
        E_y = (sin_pos_theta*sin_pos_phi*E_r)+(cos_pos_theta*sin_pos_phi*E_theta)+(cos_pos_phi*E_phi);
        E_z = (cos_pos_theta*E_r)+(-sin_pos_theta*E_theta);
    } else if (is_cylindrical){
        // Find positions in Cylindrical coordinates
        const double pos_phi = position.getphi();

        // Declare electric field in cylindrical coordinates
        const double E_rho = one_d_val;
        const double E_phi = three_d_val;

        // For efficiency, find these once and store
        const double sin_pos_phi = sin(pos_phi);
        const double cos_pos_phi = cos(pos_phi);

        // Transform to cartesian basis
        E_x = ((cos_pos_phi*E_rho)+(-sin_pos_phi*E_phi));
        E_y = ((sin_pos_phi*E_rho)+(cos_pos_phi*E_phi));
        E_z = two_d_val;
    }

    return threevector(E_x, E_y, E_z);
}

void Field_Point::find_E_field_at_point(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical, bool is_E_Field_defined){
    if (is_E_Field_defined){
        _electric_field = convert_electric_field(_one_d_val, _two_d_val, _three_d_val, _position, is_cartesian, is_spherical, is_cylindrical);
    } else {
        _electric_field = calculate_electric_field(_position, dimension, is_cartesian, is_spherical, is_cylindrical);
    }

}

threevector Field_Point::get_nearby_E_field(threevector nearby_position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical, bool is_E_Field_defined){
    threevector nearby_electric_field;
    if (is_E_Field_defined){
        nearby_electric_field = convert_electric_field(_one_d_val, _two_d_val, _three_d_val, nearby_position, is_cartesian, is_spherical, is_cylindrical);
    }
    else {
        nearby_electric_field = calculate_electric_field(nearby_position, dimension, is_cartesian, is_spherical, is_cylindrical);
    }
    return nearby_electric_field;
}

double Field_Point::get_specified_position(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical){
            double value;
            const std::string dimension_warning = "WARNING: attempting to get specified position from an invalid dimension number.";
            const std::string coord_system_warning = "WARNING: attempting to get specified position from an invalid coordinate system.";
            if (is_cartesian){
                if (dimension == 1){value = get_position_x();}
                else if (dimension == 2){value = get_position_y();}
                else if (dimension == 3){value = get_position_z();}
                else {std::cout<<dimension_warning<<std::endl;}
            } else if (is_spherical){
                if (dimension == 1){value = get_position_r();}
                else if (dimension == 2){value = get_position_theta();}
                else if (dimension == 3){value = get_position_phi();}
                else {std::cout<<dimension_warning<<std::endl;}
            } else if (is_cylindrical){
                if (dimension == 1){value = get_position_rho();}
                else if (dimension == 2){value = get_position_z();}
                else if (dimension == 3){value = get_position_phi();}
                else {std::cout<<dimension_warning<<std::endl;}
            } else {std::cout<<coord_system_warning<<std::endl;}
            return value;
}
