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
	    double value;
	    threevector position;
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

	    /** @name Mutator Functions
	     *  @brief change value of field_point
	     *  @param 
	     */
	    ///@{
	    inline void set_val(const double new_value) { value = new_value; };
	    inline void set_pos(const threevector new_pos) { position = new_pos; };
	    ///@}
	    
	    
	    /** @name Accessor functions
	     *  get the value of field point member data
	     */
	    ///@{
	    //* @var value get the field point value */
	    inline double get_value()const{ return value; };


	    //* @var position get the field point position */
	    inline threevector get_position()const{ return position; };


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

}; // end of class

#endif /*__FIELD_POINT_H_INCLUDED__ */
