/** @file Field_Point.h
 *  @brief Definition of a field_point class for use in describing Physics fields
 *  field_point Objects contain a value at a point with a three-vector (defined in file threevector.h)
 *
 *  @author Daniel Greenhouse (dg2020@ic.ac.uk (until 01.09.2020), dmgreenhouse@outlook.com (thereafter))
 *      Written immitating the style of main DiMPl author Luke Simons (ls5115@ic.ac.uk)
 *
 *  @bugs No known bugs.
 */

#ifndef __FIELD_POINT_H_INCLUDED__
#define __FIELD_POINT_H_INCLUDED__

#include "threevector.h" //!< For threevector class
#include <vector> //!< For the std::vector's, used preferentially to arrays due to easier function return
/** @class field_point
 *  @brief Definition of a field_point class for use in describing Physics fields
 *  field_point objects contain a value at a point with a three-vector (defined in file threevector.h)
 */
class Field_Point{
    protected:
	    /** @name Member Data
	     *  @brief Private data defining fields
	     *
	     * private,
	     *     value: double, value of the field at certain position
	     *     position: threevector, position in the field described as a threevector
	     */
	    ///@{
	    double _value;
	    threevector _position;
	    std::vector<double> _first_d_abc;
	    std::vector<double> _second_d_abc;
	    std::vector<double> _third_d_abc;
	    threevector _electric_field;
        threevector calculate_electric_field(threevector position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical);
	    ///@}

	    //Note to self: potentially include information about nearby points here to increase speed?

    public:
	    /** @name Constructors
	     *  @brief functions to construct a field_point class
	     */
	    ///@{
	    /** @brief Default constructor
	     *
	     * Set all member data to zero (uses default constructor of threevector)
	     */
	    Field_Point();


	    /** @brief value-position constructor
	     *  @param val defines the field point value
	     *  @param pos defines the field point position
	     */
	    Field_Point(double val, threevector pos);
        /** @brief value-position constructor
	     *  @param val defines the field point value vector
	     *  @param pos defines the field point position
	     */
	    Field_Point(threevector val, threevector pos);

	    /** @name Mutator Functions
	     *  @brief change value of field_point
	     *  @param
	     */
	    ///@{
	    inline void set_val(const double new_value) { value = new_value; };
	    inline void set_pos(const threevector new_pos) { position = new_pos; };
        inline void set_E_Field(const threevector field_values) {_electric_field = field_values};
	    ///@}

	    void add_quadratic_fit(std::vector<std::vector<double>> abc_details);
	    void find_E_field_at_point(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical);
        threevector get_nearby_E_field(threevector nearby_position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical);

	    /** @name Accessor functions
	     *  get the value of field point member data
	     */
	    ///@{
	    //* @var value get the field point value */
	    inline double get_value()const{ return value; };


	    //* @var position get the field point position */
	    inline threevector get_position()const{ return position; };

	    //* @var value the field point position at a specified dimension/coordinate system.
        // @param dimension: the dimension number being looked at
        // @param is_cartesian: if the coordinate system used is cartesian
        // @param is_spherical: if the coordinate system used is spherical
        // @param is_cylindrical: if the coordinate system used is cylindrical
        // */
	    inline double get_specified_position(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical)const{
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
                else if (dimension == 2){value = get_position_theta());}
                else if (dimension == 3){value = get_position_phi();}
                else {std::cout<<dimension_warning<<std::endl;}
            } else if (is_spherical){
                if (dimension == 1){value = get_position_rho();}
                else if (dimension == 2){value = get_position_z());}
                else if (dimension == 3){value = get_position_phi();}
                else {std::cout<<dimension_warning<<std::endl;}
            } else {std::cout<<coord_system_warning<<std::endl;}
            return value; };
	    ///@}


	    //* @var position_x get the field point x position */
	    inline double get_position_x()const{ return position.getx(); };
	    ///@}

	    //* @var position_y get the field point y position */
	    inline double get_position_y()const{ return position.gety(); };
	    ///@}

	    //* @var position_z get the field point z position */
	    inline double get_position_z()const{ return position.getz(); };
	    ///@}

	    //* @var position_theta get the field point theta position */
	    inline double get_position_theta()const{ return position.gettheta(); };
	    ///@}

	    //* @var position_phi get the field point phi position */
	    inline double get_position_phi()const{ return position.getphi(); };
	    ///@}

	    //* @var position_phi get the field point rho position */
	    inline double get_position_rho()const{ return position.getrho(); };
	    ///@}

	    //* @var position_phi get the field point r position */
	    inline double get_position_r()const{ return position.mag3(); };
	    ///@}

	    //* @var electric field at this point */
	    inline threevector get_electric_field()const{ return _electric_field; };
	    ///@}

}; // end of class

#endif /*__FIELD_POINT_H_INCLUDED__ */
