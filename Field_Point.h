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
#include "Constants.h"

using namespace dimplconsts;

/** @class field_point
 *  @brief Definition of a field_point class for use in describing Physics fields
 *  field_point objects contain a value at a point with a three-vector (defined in file threevector.h)
 */
class Field_Point{
    protected:
        ///@{
	    /** @name Private Member Data
	     *  @brief Private data defining field points
	     *
	     *     _value: double, scalar value of the field at this field point
         *             position. This is typically the potential value.
         *             This is only used when constructing according to
         *             potential values.
	     *     _one_d_val: double, component of the field value vector pointing
         *             along the first dimension at this field point position.
         *             This is typically the electric field vector component
         *             pointing along the first dimension.
         *             This is only used when constructing according to
         *             electric field values.
	     *     _two_d_val: double, component of the field value vector pointing
         *             along the second dimension at this field point position.
         *             This is typically the electric field vector component
         *             pointing along the second dimension.
         *             This is only used when constructing according to
         *             electric field values.
	     *     _three_d_val: double, component of the field value vector
         *             pointing along the third dimension at this field point
         *             position. This is typically the electric field vector
         *             component pointing along the third dimension.
         *             This is only used when constructing according to
         *             electric field values.
	     *     _position: threevector, position in the overall field at which
         *             this field point is found.
	     *     _first_d_abc: std::vector<double>, vector of three doubles
         *             containing details of the quadratic fit describing the
         *             potential fit according to neighbouring points in the
         *             first dimension. Three doubles (a,b and c respectivley)
         *             describe the potential, y, according to y = a x^2 + b x *             + c where x is the position along the first dimension.
         *             This is only used when constructing according to
         *             potential values.
	     *     _second_d_abc: std::vector<double>, vector of three doubles
         *             containing details of the quadratic fit describing the
         *             potential fit according to neighbouring points in the
         *             second dimension. Three doubles (a,b and c respectivley)
         *             describe the potential, y, according to y = a x^2 + b x *             + c where x is the position along the second dimension
         *             This is only used when constructing according to
         *             potential values.
	     *     _third_d_abc: std::vector<double>, vector of three doubles
         *             containing details of the quadratic fit describing the
         *             potential fit according to neighbouring points in the
         *             third dimension. Three doubles (a,b and c respectivley)
         *             describe the potential, y, according to y = a x^2 + b x *             + c where x is the position along the third dimension
         *             This is only used when constructing according to
         *             potential values.
	     *     _electric_field: threevector, the electric field vector at the
         *             field point position. This is either directly provided
         *             through the constructor or calculated through the
         *             potential value and abc values.
	     */
	    double _value;
        double _one_d_val;
        double _two_d_val;
        double _three_d_val;
	    threevector _position;
	    std::vector<double> _first_d_abc;
	    std::vector<double> _second_d_abc;
	    std::vector<double> _third_d_abc;
	    threevector _electric_field;
    	///@}

        ///@{
	    /** @name Member Methods
	     *  @brief Private methods utilised when defining field points
	     *
	     */
        ///@{
        /** @name calculate_electric_field
         *  @brief used to calculate the electric field
         *         according to this field point from the potential abc values
         *         at a certain nearby position. The electric field is found as
         *         the negative gradient of the potential and the different
         *         geometries accounted for.
         *  @param position: the position at which the electric field is to be
         *             calculated. This does not have to be the position of the
         *             field point itself.
         *  @param dimension: the number of dimensions being described in this
         *             field.
         *  @param is_cartesian: if this field is described in the cartesian
         *             coordinate system.
         *  @param is_spherical: if this field is described in the spherical
         *             coordinate system.
         *  @param is_cylindrical: if this field is described in the
         *             cylindrical coordinate system.
         *  @return electric field at a position according to this field point
         *          (threevector)
         */
        threevector calculate_electric_field(threevector position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical);
        ///@}
        ///@{
        /** @name convert_electric_field
         *  @brief used to convert an electric field described in a certain
         *         coordinate system to the cartesian coordinate system
         *         as used by DiMPl.
         *  @param one_d_val: the electric field component at this field point
         *             pointing along the first dimension (in whichever
         *             coordinate system the field is being described in).
         *  @param two_d_val: the electric field component at this field point
         *             pointing along the second dimension (in whichever
         *             coordinate system the field is being described in).
         *  @param three_d_val: the electric field component at this field point
         *             pointing along the third dimension (in whichever
         *             coordinate system the field is being described in).
         *  @param position: the position at which the electric field is to be
         *             calculated. This does not have to be the position of the
         *             field point itself.
         *  @param is_cartesian: if this field is described in the cartesian
         *             coordinate system.
         *  @param is_spherical: if this field is described in the spherical
         *             coordinate system.
         *  @param is_cylindrical: if this field is described in the
         *             cylindrical coordinate system.
         *  @return electric field at a position according to this field point
         *          converted to be used in a Cartesian Coordinate system
         *           (threevector)
         */
	    threevector convert_electric_field(double one_d_val, double two_d_val, double three_d_val, threevector position, bool is_cartesian, bool is_spherical, bool is_cylindrical);
        ///@}
	    ///@}

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

	    /** @brief scalar value-position constructor
         *         NOTE: this constructor requires the later use of the
         *              add_quadratic_fit method to subsequently find the
         *              electric field.
	     *  @param val: defines the field point value
	     *  @param pos: defines the field point position
	     */
	    Field_Point(const double val, threevector pos);

        /** @brief vector value-position constructor.
         *         NOTE: this constructor sets the scalar value to be 0 and so
         *              only suitable to be used when the full electric field
         *              vector is provided at construction.
	     *  @param one_d_val: defines the field point vector value along the
         *             first dimension
 	     *  @param two_d_val: defines the field point vector value along the
         *             second dimension
  	     *  @param three_d_val: defines the field point vector value along the
         *             third dimension
	     *  @param pos: defines the field point position
         *  @param is_cartesian: whether the field is described in Cartesian
         *             form or not
         *  @param is_spherical: whether the field is described in Spherical
         *             form or not
         *  @param is_cylindrical: whether the field is described in Cylindrical
         *             form or not
	     */
	    Field_Point(double one_d_val, double two_d_val, double three_d_val, threevector pos, bool is_cartesian, bool is_spherical, bool is_cylindrical);

	    /** Mutator Functions
	     *  @brief change value of field_point member data
	     */
	    ///@{
        /** @name set_val
         *  @brief sets the scalar value of the field point to a specified new
         *         value
         *  @param new_value the new scalar value assigned to the field point
         */
	    inline void set_val(const double new_value) { _value = new_value; };
        /** @name set_E_Field
         *  @brief sets the vector value of the field point to a specified new
         *         value
         *  @param field_values the new vector value assigned to the field point
         */
	    inline void set_E_Field(const threevector field_values) {_electric_field = field_values;};
	    ///@}

        /** Helper Functions
	     *  @brief public methods used by the class
	     */
	    ///@{
        /** @name add_quadratic_fit
         *  @brief adds the quadratic fit details for the scalar value
         *         description of the field point. If the electric field is
         *         being calculated from potential points then this method must
         *         be called.
         *  @param abc_details the abc details of this field point. These are
         *         contained within three vectors each containg three values.
         *         Each vector of three values have values in the order {a,b,c}
         *         (for the fit y=ax^2+bx+c) and each of the three vectors is
         *         for a separate dimension (in the order
         *         {first_dimension_vector, second_dimension_vector,
         *         third_dimension_vector}).
         */
	    void add_quadratic_fit(std::vector<std::vector<double>> abc_details);

        /** @name find_E_field_at_point
         *  @brief finds the electric field value at this field point position.
         *         If the electric field at the point is known, then it is
         *         converted such that it can correctly be interpreted in
         *         Cartesian coordinates for subsequent DiMPl usage.
         *         If the electric field at the point is not known then it is
         *         calculated but note that the abc details must have been
         *         provided.
         *         NOTE: this method finds the electric field at the point and
         *             assigns it to the member data. Consequently it must be
         *             called for the _electric_field member to be subsequently
         *             used.
         *  @param dimension the number of dimensions that the field is being
         *         described in
         *  @param is_cartesian whether the field is being described in a
         *         Cartesian coordinate system or not.
         *  @param is_spherical whether the field is being described in a
         *         Spherical coordinate system or not.
         *  @param is_cylindrical whether the field is being described in a
         *         Cylindrical coordinate system or not.
         *  @param is_E_Field_defined whether the electric field is being
         *         directly provided or not.
         */
	    void find_E_field_at_point(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical, bool is_E_Field_defined);
        ///@}

	    /** @name Accessor functions
	     *  get the value of field point member data
	     */
	    ///@{

        /** @name get_nearby_E_field
         *  @brief returns the electric field value according to this field
         *         point at a nearby position.
         *         If the electric field at the point is known, then this field
         *         is assumed not to deviate much to the nearby point and so it
         *         is converted such that it can correctly be interpreted in
         *         Cartesian coordinates for subsequent DiMPl usage.
         *         If the electric field at the point is not directly known
         *         then the field's abc details are used to find the electric
         *         field at the nearby position itself. Note that for this to be
         *         calculated the abc details must have been provided.
         *  @param nearby_position the nearby position at which the Electric
         *         field is to be found.
         *  @param dimension the number of dimensions that the field is being
         *         described in
         *  @param is_cartesian whether the field is being described in a
         *         Cartesian coordinate system or not.
         *  @param is_spherical whether the field is being described in a
         *         Spherical coordinate system or not.
         *  @param is_cylindrical whether the field is being described in a
         *         Cylindrical coordinate system or not.
         *  @param is_E_Field_defined whether the electric field is being
         *         directly provided or not.
         *  @return nearby_electric_field the electric field value at a nearby
         *         position as according to this field point. This field has
         *         been converted to be correctly interpreted by Cartesian DiMPl
         */
        threevector get_nearby_E_field(threevector nearby_position, int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical, bool is_E_Field_defined);

	    /** @name get_value
         *  @brief returns the scalar value of the field point.
         *  @return _value the field point's scalar value
         */
	    inline double get_value()const{ return _value; };

	    /** @name get_specified_position
         *  @brief returns the position value of the field in the dimension
         *        specified. This allows for all coordinate systems and follows
         *        the form:
         *        Cartesian: (x,y,z) for (1st,2nd and 3rd dimensions)
         *        Spherical: (r, theta , phi) for (1st,2nd and 3rd dimensions)
         *        Cylindrical: (rho, z , phi) for (1st,2nd and 3rd dimensions)
         *  @param dimension: the dimension in which the position value is to
         *         be Found
         *  @param is_cartesian: whether the field is described in Cartesian
         *         form or not
         *  @param is_spherical: whether the field is described in Spherical
         *         form or not
         *  @param is_cylindrical: whether the field is described in Cylindrical
         *         form or not
         * @return value: the position value at the specified dimension
         *         according to the specified coordinate system.
         */
	    double get_specified_position(int dimension, bool is_cartesian, bool is_spherical, bool is_cylindrical);

        /** @name get_position_x
         *  @brief uses the field point's member data _position to get the
         *        field point's x position value
         *  @return the field point's x position value
         */
	    inline double get_position_x()const{ return _position.getx(); };

        /** @name get_position_y
         *  @brief uses the field point's member data _position to get the
         *        field point's y position value
         *  @return the field point's y position value
         */
	    inline double get_position_y()const{ return _position.gety(); };

        /** @name get_position_z
         *  @brief uses the field point's member data _position to get the
         *        field point's z position value
         *  @return the field point's z position value
         */
	    inline double get_position_z()const{ return _position.getz(); };

        /** @name get_position_theta
         *  @brief uses the field point's member data _position to get the
         *        field point's theta position value
         *  @return the field point's theta position value
         */
	    inline double get_position_theta()const{ return _position.gettheta(); };

        /** @name get_position_phi
         *  @brief uses the field point's member data _position to get the
         *        field point's phi position value
         *  @return the field point's phi position value
         */
	    inline double get_position_phi()const{ return _position.getphi(); };

        /** @name get_position_rho
         *  @brief uses the field point's member data _position to get the
         *        field point's rho position value
         *  @return the field point's rho position value
         */
	    inline double get_position_rho()const{ return _position.getrho(); };

        /** @name get_position_r
         *  @brief uses the field point's member data _position to get the
         *        field point's r position value
         *  @return the field point's r position value
         */
	    inline double get_position_r()const{ return _position.mag3(); };

        /** @name get_electric_field_at_point
         *  @brief uses the field point's member data _electric_field to get the
         *        field point's electric field vector. Note this must have been
         *        found first using the find_E_field_at_point method.
         *  @return the field point's electric field value at the point in
         *        vector form such that DiMPl can interpret it in Cartesian form
         */
	    inline threevector get_electric_field_at_point()const{ return _electric_field; };
	    ///@}

}; // end of class

#endif /*__FIELD_POINT_H_INCLUDED__ */
