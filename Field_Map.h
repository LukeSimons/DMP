/** @file Field_Map.h
 *  @brief Definition of a Field_Map class for use in describing Physics fields
 *  Field_Map Objects contain details of a field comprised of many Field_Point objects
 *
 *  @author Daniel Greenhouse (dg2020@ic.ac.uk (until 01.09.2020), dmgreenhouse@outlook.com (thereafter))
 *      Written immitating the style of main DiMPl author Luke Simons (ls5115@ic.ac.uk)
 *
 *  @bugs No known bugs.
 */

#ifndef __FIELD_MAP_H_INCLUDED__
#define __FIELD_MAP_H_INCLUDED__

#include "threevector.h" //!< For threevector class
#include "Field_Point.h" //!< For Field_Point class

#include <iostream>
#include <algorithm>
#include <cctype>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <math.h> //!< For log 2

/** @class Field_Map
 *  @brief Definition of a Field_Map class for use in describing Physics fields
 *  Field_Map Objects contain details of a field comprised of many Field_Point objects
 *
 */
class Field_Map{
    protected:
        /** @name Member Data
	 *  @brief Private data defining field maps
	 */
	std::vector<std::string> _line_vector; // used to store contents of each line in text file, vector used so that the text file can be an arbitrary size. Vector populated from constructor
	int _num_data_lines; // the number of lines containing data of the field in the text file. Value found in the constructor
	int _num_data_points_per_line; // the number of data points at each field point (eg. for 3D field there will be four points (fourth being the value of the field at that point);

	// Header Details:
	//     details necessary header structure for the text file
	const std::string _summary_delim = "Summary:\t";
	const std::string _dimension_delim = "Number of Dimensions:\t";
	const std::string _coord_delim = "Coordinate System:\t";
	const std::string _ordered_delim = "Ordered:\t";
        const int _num_lines_in_header = 4; // the number of lines in the header of the text file

	std::string _summary_line; // found in the constructor
	int _dimension_val; // the number of dimensions used to describe the electric field. Found in the constructor
	bool _ordered_truth; // found in the constructor
        std::string _coord_type; // Coordinate system being used for electric field description. Found in the constructor. See below for details:
            /** "Cart_1" = 1D Cartesian,
             *  "Cart_2" = 2D Cartesian,
             *  "Cart_3" = 3D Cartesian,
             *  "Polar_1" = 1D Polar (Radial direction only),
             *  "Polar_2" = 2D Polar,
             *  "Spherical_3" = 3D Spherical,
             *  "Cylindrical_3" = 3D Cylindrical,
             *  Note to self: this could be better written as an enumeration
             */
        // Data Details:
	//     details necessary data structure for the text file
	const std::string _data_delim = "\t"; // the delimiter between data values
        void convert_data_line_to_multi_D_array();
        void text_to_string_vector(std::string field_file_name);
	std::string split_after_delim(std::string delim, std::string str_to_split);
	void find_header_details();
	bool string_to_bool(std::string string);
	int _total_one_d_count;
	int _total_two_d_count;
	int _total_three_d_count;
	std::vector<std::vector<std::vector<Field_Point>>> _Field_Map_Final;
	void partition_three_D_grid();
	std::vector<int> partition(int start, int end);
	std::vector<int> find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder);
	void find_per_dimension(int dimension);
	std::vector<std::vector<int>> _one_d_pos_holder;
	std::vector<std::vector<double>> _one_d_value_holder;
	std::vector<std::vector<int>> _two_d_pos_holder;
	std::vector<std::vector<double>> _two_d_value_holder;
	std::vector<std::vector<int>> _three_d_pos_holder;
	std::vector<std::vector<double>> _three_d_value_holder;

    public:
	/** @name Constructors
	 *  @brief constructs a Field_Map class using a specfied string detailing the location of a custom Field Map
	 */
	///@{
	/** @brief Default Constructor
	 *  @param field_file_name std::string detailing the location of a custom Field Map
	 */
	Field_Map(std::string field_file_name);

	double find_approx_value(threevector point);
}; //end of class

#endif /*__FIELD_MAP_H_INCLUDED__ */
