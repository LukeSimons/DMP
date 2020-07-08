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

