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
#include <ctime> //!< For timer, Note to self: remove this
#include <cstdio> //!< For clock time conversion, Note to self: remove this

#include "Constants.h"

using namespace dimplconsts;

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

        ///@{
    	/** @brief Member Variables
        *          Class wide member variables are denoted by an underscore prefix.
        *          Many member variables are found during the constructor processes.
        *   @param _DEBUG: bool; used to optionally print additional information if required.
        *   @param _line_vector: std::vector<std::string>; used to store contents of each line in text file, vector used so that the text file can be an arbitrary size. Vector populated from constructor
        *   @param _Field_Map_Final: std::vector<std::vector<std::vector<Field_Point>>>; used to store the contents of each data point as a Field_Point in a 3D grid. std::vector is preferred to arrays since they are more easily returned by methods.
        *   @param _one_d_pos_holder: std::vector<std::vector<int>>; Contains information of the partitioning in the first dimension. The first line contains a vector of three values: the first, middle and last position. Each subsequent line contains a vector containing additional positions to the line above where these additional positions are the midpoints between each position in the line above. The final line is every single position in the _Field_Map_Final for that dimension.
        *   @param _one_d_value_holder: std::vector<std::vector<double>>; the value of the position in that dimension corresponding to the index in _one_d_pos_holder
        *   @param _two_d_pos_holder: as with _one_d_pos_holder but for the second dimension
        *   @param _two_d_value_holder: as with _one_d_value_holder but for the second dimension
        *   @param _three_d_pos_holder: as with _one_d_pos_holder but for the third dimension
        *   @param _three_d_value_holder: as with _one_d_value_holder but for the third dimension
        */
        const bool _DEBUG = true;
        const bool _CALCULATE_NEARBY_E_FIELD_ON_CURVES = true;
	    std::vector<std::string> _line_vector;
    	std::vector<std::vector<std::vector<Field_Point>>> _Field_Map_Final;
    	std::vector<std::vector<int>> _one_d_pos_holder;
    	std::vector<std::vector<double>> _one_d_value_holder;
    	std::vector<std::vector<int>> _two_d_pos_holder;
    	std::vector<std::vector<double>> _two_d_value_holder;
    	std::vector<std::vector<int>> _three_d_pos_holder;
    	std::vector<std::vector<double>> _three_d_value_holder;
	std::vector<std::vector<int>> _num_times_outside_domain{3, std::vector<int> (2, 0)};
        std::vector<std::vector<double>> _two_D_array;
        ///@{
        /** Header Details:
        *       Necessary Header structure in the text file
        *   @param _summary_delim: std::string; the prefix to the Summary information in the text file header
        *   @param _dimension_delim: std::string; the prefix to the Dimension Number information in the text file header
        *   @param _coord_delim: std::string; the prefix to the Coordinate Type information in the text file header
        *   @param _ordered_delim: std::string; the prefix to the Ordered truth information in the text file header
        *   @param _num_lines_in_header: int; the number of lines in the header of the text file
        *   @param _summary_line: std::string; holds the information of the Summary line as a single string
        *   @param _dimension_val: int; the number of dimensions used to describe the field
        *   @param _ordered_truth: bool; whether the provided text file is in the correctly ordered form or not.
        *   @param _coord_type: std::string; the coordinate system being used for field description. See below for details:
        *         * @param _CART_1: Cartesian_1" = 1D Cartesian (extended to 3D Cartesian),
        *         * @param _CART_2: "Cartesian_2" = 2D Cartesian (extended to 3D Cartesian),
        *         * @param _CART_3: "Cartesian_3" = 3D Cartesian,
        *         * @param _SPHER_1: "Spherical_1" = 1D Radial (extended to Spherical),
        *         * @param _SPHER_2: "Spherical_2" = 2D Polar (extended to Spherical),
        *         * @param _SPHER_3: "Spherical_3" = 3D Spherical,
        *         * @param _CYL_1: "Cylindrical_1" = 1D Radial (extended to Cylindrical),
        *         * @param _CYL_2: "Cylindrical_2" = 2D Polar (extended to Cylindrical),
        *         * @param _CYL_3: "Cylindrical_3" = 3D Cylindrical,
        *         *  Note to self: this could be better written as an enumeration
        *   @param _is_cartesian: boolean; the 3D coordinate system being used for the field description. If only 1 or 2 dimensions provided, the additional dimensions will be flooded with "0.0"s in Cartesian fashion.
        *       Cartesian systems are in the form: (x,y,z)
        *   @param _is_spherical: boolean; the 3D coordinate system being used for the field description. If only 1 or 2 dimensions provided, the additional dimensions will be flooded with "0.0"s in Spherical fashion.
        *       Spherical systems are in the form: (r,theta,phi). Note this change in convention is since typically phi will be symmetrical
        *   @param _is_cylindrical: boolean; the 3D coordinate system being used for the field description. If only 1 or 2 dimensions provided, the additional dimensions will be flooded with "0.0"s in Cylindrical fashion.
        *       Cylindrical systems are in the form: (rho,z,phi). Note this change in convention is since typically phi will be symmetrical
        */
        const std::string _summary_delim = "Summary:\t";
	    const std::string _dimension_delim = "Number of Dimensions:\t";
	    const std::string _coord_delim = "Coordinate System:\t";
	    const std::string _ordered_delim = "Ordered:\t";
	    const std::string _is_E_Field_defined_delim = "Is Electric Field Defined:\t";
	    const std::string _end_delim = ";";
        const int _num_lines_in_header = 5;
        std::string _summary_line;
	    int _dimension_val;
	    bool _ordered_truth;
	    bool _is_E_Field_defined_truth;
        std::string _coord_type;
        const std::string _CART_1 = "Cartesian_1";
        const std::string _CART_2 = "Cartesian_2";
        const std::string _CART_3 = "Cartesian_3";
        const std::string _SPHER_1 = "Spherical_1";
        const std::string _SPHER_2 = "Spherical_2";
        const std::string _SPHER_3 = "Spherical_3";
        const std::string _CYL_1 = "Cylindrical_1";
        const std::string _CYL_2 = "Cylindrical_2";
        const std::string _CYL_3 = "Cylindrical_3";
        bool _is_cartesian = false;
        bool _is_spherical = false;
        bool _is_cylindrical = false;
        ///@}
        ///@{
        /** Data Details:
        *       Necessary Data structure in the text file
        *   @param _data_delim: std::string; delimiter between data values
        *   @param _num_data_lines: int; the number of lines containing data of the field in the text file. Value found in the constructor
        *   @param _num_data_points_per_line: int; the number of data points at each field point (eg. for 3D field there will be four points (fourth being the value of the field at that point);
        *   @param _total_one_d_count: int; the total number of data points used to describe the field in the first dimension
        *   @param _total_two_d_count: int; the total number of data points used to describe the field in the second dimension
        *   @param _total_three_d_count: int; the total number of data points used to describe the field in the third dimension
        */
	    const std::string _data_delim = "\t";
	    int _num_data_lines;
	    int _num_data_points_per_line;
    	int _total_one_d_count;
    	int _total_two_d_count;
    	int _total_three_d_count;
        ///@}
        ///@}

    	/** @name Member Methods
        *   @brief Class wide member methods to be used within the class.
        *          Many member methods are called and managed during the constructor processes.
	*/
        ///@{
        /** @brief text_to_string_vector:
         *      Converts a text file to a vector containing contents of each line.
         *      Each element of the vector is a line of the text file in otder.
         *      The Vector is of initially arbitrary size and scales to the size of the text file.
         *  @param field_file_name : the name of the text file
         *  @param _line_vector: method wide container of the vector of strings.
         */
         void text_to_string_vector(std::string field_file_name);

         /** @brief find_header_details:
          *      Takes the header information (stored as strings) and extracts the header details from them
          *  @param _line_vector: method wide container of the vector of strings.
          *  @submethod split_after_delim: used to extract useful information from string
          *  @submethod string_to_bool: used to convert the string to a boolean
          */
     	 void find_header_details();

         /** @brief split_after_delim:
          *      Takes a string and returns the substring after a provided delimiter
          *  @param delim: std::string; a string denoting what is expected at the start of the main string.
          *  @param str_to_split: std::string; the main string to be split. This main string shall start with the delimiter and end with the important information.
          *  @return the string section containing the important information
          * Note to self: this can be improved by providing an option to throw an error if the delimiter is not found.
          */
	     std::string split_after_delim(std::string delim, std::string str_to_split);

         /** @brief string_to_bool:
          *      Takes a string and returns a boolean
          *  @return the corresponding boolean
          */
	     bool string_to_bool(std::string string);

         /** @brief convert_data_line_to_multi_D_array:
          *      Take raw data (in an ordered fashion) and transform it to a multi dimensional array.
          *  Procedure:
          *      Begins by recording the data dimensions.
          *      Takes the vector of strings and transforms it to a 2D array of the data of shape [I][J]. I: Value, 1D position, 2D position (optional), 3D position (optional). J: next set at next ordered position.
          *      Next, for the first dimension runs through the 2D array line by line and counts whilst the 1D position number is ascending. e.g. if the 1D position numbers go {-1,0,1,2,3,4,-1,0,1,2,3,4,-1...} then the count will be 6.
          *      Optionally (depending on if there are more than 1 dimensions specified in the text file), repeat the procedure for the 1D position number count but with the 2D position. Since each increment in the 2D position occurs after a full count of the 1D position, these are checked after i*one_d_count.
          *      Optionally (depending on if there are more than 2 dimensions specified in the text file), repeat the procedure for the 1D position number count but with the 3D position. Since each increment in the 3D position occurs after a full count of both the 1D position and 2D position, these are checked after i*one_d_count*two_d_count.
          *      Populate a 3 dimensional grid with the values, where the grid spans real space. Each point is a Field_Point object. Conversion from the 2D array to the 3D grid is performed using the formula: (k*one_d_count*two_d_count) + j * one_d_count + i where k is the third dimension position, j is the second dimension position and i is the first dimension position. The Field_Point objects recognise a vector as being either cartesian, spherical or cylindrical depending on previously found type specified in file header.
          *  @param two_D_array: double; contains the data of the text file, without altering the structure of it. If only 1 or 2 dimensions are defined, the additional (up until 3) are flooded with 0's. Eg. if 1 dimension is defined across 20 values, there shall still only be 20 data points, but now these contain a position of 0 and 0 in the other two dimensions.
          */
         void convert_data_line_to_multi_D_array();

         /** @brief partition_three_D_grid:
          *      Used to handle the partitioning of the 3D grid of Field_Point objects. Works according to the number of dimensions being used to describe the field - it only partitions those dimensions that are described.
          *      If the dimension is not described, the `partition' is simply the single point.
          *  @submethod: find_per_dimension; actually does the partitioning
          */
	     void partition_three_D_grid();

         /** @brief find_per_dimension:
          *      Used to partition a Grid along a single dimension.
          *  Procedure:
          *      Specifies which dimension is being partitioned and manages that dimension of the 3D grid accordingly.
          *      Uses a `pos_holder' to account for the index of the grid point (the value 3 in a grid 2,3,4,5 shall have an index/pos of 1) and a `val_holder' to account for the corresponding value of interest (the position/distance in that dimension at the correspondingpos/index).
          *      The grid in that dimension is split up into sections (@param num_secs) where each section contains an additional level of granularity in that dimension.
          *      The first line is the two end points and the central points (0, n/2, n). For the pos, if n/2 is not an integer then it is rounded down. This is managed by the @submethod partition.
          *      The second line is the above line but with additional midpoints between all points (0, n/4, n/2, 3n/4, n). Again, non-integer pos's are rounded down.
          *      The third to penultimate line follow the same procedure (adding the midpoints each time)
          *      The final line is every single point.
          *      The corresponding class-wide member variable position and value holders for that dimension are updated accordingly
          *  @param dimension: int; the dimension in which the values are being partitioned
          */
     	 void find_per_dimension(int dimension);

         /** @brief partition:
          *      Used to find the midpoint between two positions.
          *  @param start: int; the index of the first position
          *  @param end: int; the index of the final position
          *  @param midway: int; the index of the midway position between the start and end. This will always return the rounded-down integer.
          *  @return vector containing the start, midaway and end index position
          */
	     std::vector<int> partition(int start, int end);

         /** @brief find_val_from_field_point:
          *      Used to get the corresponding value from a field point.
          *      The ordering of the cartesian, spherical and cylindrical are (x,y,z), (r, theta, phi) and (rho, z, phi) for dimension 1,2,3 respectively.
          *  @param index: int; the index of the Field_Point being looked at
          *  @param dimension: int; the dimension being looked at
          *  @return value: double; the found corresponding value
          */
         double find_val_from_field_point(int index, int dimension);
         /** @brief add_quadratic_fit_details_to_all_points:
          *      Used to manage the finding of quadratic fits between points for all points in the Field_Map.
          *      Also calls a method on each Field_Point to record the fit details associated with the point.
          *      Furthermore calls a method on each Field_Point to find the electric field at that point.
          *  @param pos: int[3]; the index of the position in the grid in terms of {first index position, second index position, third index position}
          *  @param fits: std::vector<std::vector<double>>; each line refers to one of the three dimensions and along each line is the a,b and c values referring to the fit v = a*s^2+b*s+c where v is the value of the field at that point (a potential) and s is the position in space along that dimension.
          *  @submethod find_fits
          */
 	     void add_quadratic_fit_details_to_all_points();

         /** @brief find_fits:
          *    Specifies the neighbouring points in each dimension required to find the fit.
          *    Procedure
          *        Checks if the point being looked at is an edge point in which case it is treated differently (the neighbouring point is set as the point itself).
          *        Finds the neighbouring points (1 up and 1 down) in each of the three dimensions
          *        Retrieves the corresponding value at these points and the position in that dimension
          *        Calls the @submethod calc_quad_fit to find the relationship of the central point to the neighbouring points
          *    @return abc_s: std::vector<std::vector<double>>; the abc values for the quadratic fit for each of the three dimensions (stored as each line)
          */
     	 std::vector<std::vector<double>> find_fits(int position[3]);

         /** @brief calc_quad_fit:
          *    Calculates the fit equation corresponding to a point according to its relationship with neighbouring points. The fit is of the form y=ax^2+bx+c (quadratic)
          *    Procedure
          *        Creates a set of three linear equations specifying the coordinates that a quadratic curve must fit through and solves them using a 3x3 matrix
          *        Checks if the point being looked at is an edge point in which the linear equations shall have repeated points, hence these are solved accordingly.
          *    @param p1: double[2]; the position (in the dimension of interest) and value of the low point.
          *    @param p2: double[2]; the position (in the dimension of interest) and value of the mid point.
          *    @param p3: double[2]; the position (in the dimension of interest) and value of the high point.
          *    @param _is_repeated_low_point: bool; if the mid point is on the lower edge position of that dimension, then the low point shall just be a repitition of this
          *    @param _is_repeated_high_point: bool; if the mid point is on the higher edge position of that dimension, then the high point shall just be a repitition of this
          *    @return Y; std::vector<double>; the details of the fit written as {a,b,c}
          */
         std::vector<double> calc_quad_fit(double p1[2], double p2[2], double p3[2]);

         /** @brief two_by_two_mat_det:
          *    Calculates the determinant of a 2x2 matrix
          *    @param upper_left: double; the 00 value of the matrix
          *    @param upper_right: double; the 01 value of the matrix
          *    @param bottom_left: double; the 10 value of the matrix
          *    @param bottom_right: double; the 11 value of the matrix
          *    @return the determinant of the matrix (double)
          */
         double two_by_two_mat_det(double upper_left, double upper_right, double bottom_left, double bottom_right);

         /** @brief find_closest:
          *    Uses the partitioned matrices of the index positions and values (the position in a certain dimension and coordinate system) to find the index of where the position being sought is.
          *    Procedure:
          *        Begin by looking at the values in each line.
          For the first line, if the true value is above the midpoint then on the line beneath the next midpoint should be considered.
          This is enacted by keeping track of a @param`test_point' which is shifted accordingly for each subsequent line.
          This is continued up until the line before the penultimate line (since the final line is spaced differently)
          *        If the value is found to be below the `test_point' in the penultimate line, then the index below the `test_point'  is the base and the `test_point' is the upper. Conversely, if the value is found to be above the `test_point' in the penultimate line then the `test_point' is the base and the index above the `test_point'  is the upper.
          *        Next the spacing between the upper and base positions are found (@params num_to_check) and each of these position indices are checked on the final line to check which is the closest position and which is the second closest.
          *    @param value: double; the position being sought.
          *    @param pos_holder: std::vector<std::vector<int>>; the partitioned position holder matrix
          *    @param value_holder: std::vector<std::vector<double>>; the partitioned value holder matrix
          *    @return closest_positions: std::vector<int> the indices of the closest and second closest positions.
          */
	     std::vector<int> find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder, int dimension);

         void check_coordinate_domain_symmetry(int dimension_to_check);
         void check_coordinate_domain_completeness(int dimension_to_check);
         std::vector<double> find_complicit_domain_limits(bool is_spherical_injection, bool is_cylindrical_injection, double z_min_dimpl_pars, double z_max_dimpl_pars, double max_radius_dimpl_pars);
         void check_dust_grain_map_coverage(double a1, double a2, double a3);
         std::string get_coordinate_name(int dimension);
         void check_domain(bool is_spherical_injection, bool is_cylindrical_injection);

    public:
	/** @name Constructors
	 *  @brief constructs a Field_Map class using a specfied string detailing the location of a custom Field Map
	 */
	///@{
	/** @brief Default Constructor. Manages entire setup of a Field_Map object through the calling of various methods each of which is essential for the full description of a Field_Map.
     *          Procedure:
     *            Call method to transform text file to vector of strings
     *            Call method to find the header details for the Field Map construction
     *            Call method to transform vector of strings to multi dimensional array of Field Point objects
     *            Call method to partition the multi dimensional array for easy position searching
     *            Calls method to find the relationship between neighbouring Field_Points to create a quadratic fit between them.
     *   @param field_file_name: std::string; the name of the file containing the field details.
     */
	Field_Map(std::string field_file_name);
    ///@}

    /** @name Public methods
    */
    ///@{
    /** @brief find_approx_value:
    *        Searches in the Field_Map for the approximate field value at the position provided. The position provided does not need to be an exact point described in the Field_Map.
    *   @param point: threevector; the position threevector of the point at which the field value is required.
    *   @return threevector; approximate electric field value at the point provided
    */
	threevector find_approx_value(threevector point, bool is_print);
    ///@}
	inline std::vector<std::vector<int>> get_num_times_outside_domain(){return _num_times_outside_domain;};

        void check_spherical_coordinate_range(double ImpactParameter, double maxfactor);

        void check_cylindrical_coordinate_range(double lower_z_lim, double upper_z_lim);
}; //end of class

#endif /*__FIELD_MAP_H_INCLUDED__ */
