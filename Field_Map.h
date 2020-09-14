/** @file Field_Map.h
 *  @brief Definition of a Field_Map class for use in describing Physics fields
 *  Field_Map Objects contain details of a field comprised of many Field_Point objects
 *
 *  @author Daniel Greenhouse (dg2020@ic.ac.uk (until 01.09.2020), dmgreenhouse@outlook.com (thereafter))
 *      Written immitating the style of main DiMPl author Luke Simons (ls5115@ic.ac.uk)
 *
 *  @bug shrink_map_to_fit is not working, results in a segmentation fault
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
#include <sstream> //!< For string to double conversion (istringstream)

#include <fstream>

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

    	/** @name Member Variables
         *  @brief Class wide member variables are denoted by an underscore
         *      prefix. Many member variables are found during the pseudo
         *      constructor process.
         *
         *  @var _DEBUG: used to optionally print additional information if
         *      required. Many feautures of DEBUG have been removed, but has
         *      been found to be useful in conjunction with a timer to find the
         *      sticking points if the programme is running slow.
         *  @var _CALCULATE_NEARBY_E_FIELD_ON_CURVES: if true, then the
         *      programme shall find the approximate electric field according
         *      to each field point on the quadratic curve at a nearby point
         *      but if false then it shall simply approximate the nearby
         *      electric field to be the same as that at the field point itself.
         *      Runtime vs accuracy considerations should be made when choosing
         *      this toggle.
         *  @var _one_d_pos_holder: Contains information of the partitioning in
         *      the first dimension. Each element is an index position in the
         *      map rather than the true position value. The first line
         *      contains a vector of three values: the first, middle and last
         *      position. Each subsequent line contains a vector containing
         *      additional positions to the line above where these additional
         *      positions are the midpoints between each position in the line
         *      above. The final line is every single position in the
         *      _Field_Map_Final for that dimension.
         *  @var _one_d_value_holder: as with _one_d_pos_holder (see above) but
         *      with the actual position value stored in each element rather
         *      than the index of the position. (eg. if the first radial
         *      position at which there is a data point in a spherical map is
         *      0.1 then the first element of the value_holder shall be 0.1,
         *      but the first element of the pos_holder shall be 0).
         *  @var _two_d_pos_holder: as with _one_d_pos_holder (see above) but
         *      for the second dimension.
         *  @var _two_d_value_holder: as with _one_d_pos_holder (see above) but
         *      for the second dimension.
         *  @var _three_d_pos_holder: as with _one_d_pos_holder (see above) but
         *      for the third dimension.
         *  @var _three_d_value_holder:as with _one_d_pos_holder (see above) but
         *      for the third dimension.
         *  @var _line_vector: used to store contents of each line in the
         *      field map text file as a vector (each element is a line of
         *      text). vector is used so that the text file can be an
         *      arbitrary size. Vector populated from pseudo constructor.
         *  @var _Field_Map_Final: used to store the contents of each data
         *      point as a Field_Point in a 3D grid. By this vector is
         *      preferred to arrays since they are more easily returned by
         *      methods. If vectors are found to be slow for very large files,
         *      then may be worth changing this.
         *  @var _num_times_outside_domain: stores details of how many times a
         *      position was sought that was not available within the domain of
         *      the custom field map. Matrix form, each row is a new dimension
         *      (1st, 2nd, 3rd) and each column is below the domain limit (1st)
         *      or above the domain limit (2nd).
         *  @var _two_D_array: holds all of the data (not the header details)
         *      from the text file as a large two dimensional array where each
         *      row is a data point and the columns within that row are the
         *      details of that data point (potential, first dimensional
         *      position etc.). This replicates the same form as the text file.
         *  @var _eTemp: the electron temperature (in eV)
         *  @var _debye_length the Debye Length (in m) as found from the main
         *      DiMPl pars
         *  @var _dust_radius: the dust radius (in m)
         *  @var _dust_potential: the dust potential in normalised units ( phi
         *      = (-eV)/(kBTe))
         *  @var _is_spherical_injection: whether spherical injection is to be
         *      used by DimPl or not.
         *  @var _is_cylindrical_injection: whether cylindrical injection is
         *      to be used by DimPl or not.
         *  @var _a1: the dust radius 1/(semix^2) paramater used by DiMPl (in
         *      units of the dust radius)
         *  @var _a2: the dust radius 1/(semiy^2) paramater used by DiMPl (in
         *      units of the dust radius)
         *  @var _a3: the dust radius 1/(semiz^2) paramater used by DiMPl (in
         *      units of the dust radius)
         */
        const bool _DEBUG = false;
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
        double _eTemp;
        double _debye_length;
        double _dust_radius;
        double _dust_potential;
        bool _is_spherical_injection;
        bool _is_cylindrical_injection;
        double _a1;
        double _a2;
        double _a3;

        /** @name Header Details:
         *  @brief Contains information on the necessary Header structure in
         *      the text file.
         *      These should be added to and/or altered if additional
         *      header details are to be provided in future text files.
         *  @var _summary_delim: the prefix to the Summary information in the
         *         text file header
         *  @var _dimension_delim: the prefix to the Dimension Number
         *         information in the text file header
         *  @var _coord_delim: the prefix to the Coordinate Type information in
         *         the text file header
         *  @var _ordered_delim: the prefix to the Ordered truth information in
         *         the text file header
         *  @var _is_E_Field_defined_delim: the prefix to the is the electric
         *         field defined in the text file truth information in the text
         *         file header
         *  @var _is_pos_dust_radius_normalised_delim: the prefix to the is
         *        the position information in the text file normalised to the
         *        dust radius truth information in the text file header
         *  @var _is_pos_debye_length_normalised_delim: the prefix to the is
         *        the position information in the text file normalised to the
         *        debye length truth information in the text file header
         *  @var _is_potential_normalised_delim: the prefix to the is
         *        the potential information in the text file normalised to the
         *        usual scheme ((-eV)/(kBTe)) truth information in the text
         *        file header
         *  @var _is_initial_dust_pot_delim: the prefix to the is the first
         *        potential value the dust potential truth information in the
         *        text file header
         *  @var _is_initial_dust_rad_delim: the prefix to the is the first
         *        position value the dust radius truth information in the
         *        text file header
         *  @var _end_delim: the delimiter used to signal the end of a line in
         *        the header details of a text file.
         *  @var _num_lines_in_header: the number of lines in the header of the
         *         text file
         *  @var _summary_line: holds the information of the Summary line as a
         *         single string
         *  @var _dimension_val: the number of dimensions used to describe the
         *         field (1, 2 or 3)
         *  @var _ordered_truth: whether the provided text file is in the
         *        correctly ordered form or not.
         *  @var _is_E_Field_defined_truth: whether the electric field is
         *        defined in the text file or not. If it is not then it is
         *        assumed that the potential is instead defined.
         *  @var _is_pos_dust_radius_normalised_truth: whether the position
         *        values in the text file are normalised to the dust radius or
         *        not.
         *  @var _is_pos_debye_length_normalised_truth: whether the position
         *        values in the text file are normalised to the debye length or
         *        not.
         *  @var _is_potential_normalised_truth: whether the potential values
         *        in the text file are normalised to the (-eV)/(kBTe) scheme or
         *        not
         *  @var _is_initial_dust_pot_truth: whether the first potential value
         *        in the text file is the dust potential or not.
         *  @var _is_initial_dust_rad_truth: whether the first position value
         *        in the text file is the dust radius or not.
         */
         const std::string _summary_delim = "Summary:\t";
         const std::string _dimension_delim = "Number of Dimensions:\t";
         const std::string _coord_delim = "Coordinate System:\t";
         const std::string _ordered_delim = "Ordered:\t";
         const std::string _is_E_Field_defined_delim = "Is Electric Field Defined:\t";
         const std::string _is_pos_dust_radius_normalised_delim = "Positions Normalised to Dust Radius:\t";
         const std::string _is_pos_debye_length_normalised_delim = "Positions Normalised to Debye Length:\t";
         const std::string _is_potential_normalised_delim = "Potential/Electric Field Normalised via phi = (-eV)/(kBTe):\t";
         const std::string _is_initial_dust_pot_delim = "Find Dust Grain Potential from this custom field map:\t";
         const std::string _is_initial_dust_rad_delim = "Find Dust Grain Radius from this custom field map:\t";
         const std::string _end_delim = ";";
         const int _num_lines_in_header = 11;
         std::string _summary_line;
 	     int _dimension_val;
 	     bool _ordered_truth;
         bool _is_E_Field_defined_truth;
         bool _is_pos_dust_radius_normalised_truth;
         bool _is_pos_debye_length_normalised_truth;
         bool _is_potential_normalised_truth;
         bool _is_initial_dust_pot_truth;
         bool _is_initial_dust_rad_truth;

        /**
         *  @name Coordinate Sytem details:
         *  @brief Contains information on  the coordinate system being used by
         *      this map. This is a subcategory of header details since it is
         *      also extracted from the header detail section of the text file.
         *      Cartesian systems are in the form: (x,y,z)
         *      Spherical systems are in the form: (r,theta,phi)
         *      Cylindrical systems are in the form: (rho,z,phi). Note this
         *      change in convention is since phi will typically be
         *      symmetrical more often than z
         *  @var _coord_type: the coordinate system being used for field
         *       description.
         *  @var _CART_1: "Cartesian_1" = 1D Cartesian coordinate system
         *  @var _CART_2: "Cartesian_2" = 2D Cartesian coordinate system
         *  @var _CART_3: "Cartesian_3" = 3D Cartesian coordinate system
         *  @var _SPHER_1: "Spherical_1" = 1D Radial coordinate system
         *       (extended to Spherical),
         *  @var _SPHER_2: "Spherical_2" = 2D Polar coordinate system
         *       (extended to Spherical),
         *  @var _SPHER_3: "Spherical_3" = 3D Spherical coordinate *       system
         *  @var _CYL_1: "Cylindrical_1" = 1D Radial coordinate system
         *       (extended to Cylindrical),
         *  @var _CYL_2: "Cylindrical_2" = 2D Polar coordinate system
         *       (extended to Cylindrical),
         *  @var _CYL_3: "Cylindrical_3" = 3D Cylindrical coordinate
         *       system
         *  @var _is_cartesian: the 3D coordinate system being used for
         *       the field description. If only 1 or 2 dimensions
         *       provided, the additional dimensions will be flooded
         *       with "0.0"s in Cartesian fashion.
         *  @var _is_spherical: the 3D coordinate system being used for
         *       the field description. If only 1 or 2 dimensions
         *       provided, the additional dimensions will be flooded
         *       with "0.0"s in Spherical fashion.
         *  @var _is_cylindrical: the 3D coordinate system being used for
         *       the field description. If only 1 or 2 dimensions
         *       provided, the additional dimensions will be flooded
         *       with "0.0"s in Cylindrical fashion.
         */
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

        /** @name Domain Details:
         *  @brief Contains information on the domain extent of this map
         *  @var _map_z_min: the minimum extent of the map in the z-direction *      (in units of dust radius)
         *  @var _map_z_max: the maximum extent of the map in the z-direction *      (in units of dust radius)
         *  @var _map_rho_max: the maximum extent of the map in the rho-
         *      direction (radial in the x-y plane) (in units of dust radius)
         */
        double _map_z_min;
        double _map_z_max;
        double _map_rho_max;

        /** @name Data Details:
         *  @brief Contains information on the necessary Data structure in the
         *      text file
         *  @var _data_delim: delimiter between data values in the text file
         *  @var _num_data_lines: the number of lines containing data of the
         *       field in the text file. Value found during the pseudo
         *       constructor
         *  @var _num_data_points_per_line: the number of data points at each
         *        field point (eg. for 3D field using potential values rather
         *        than a defined electric field there will be four points)
         *  @var _total_one_d_count: the total number of data points used to
         *        describe the field in the first dimension
         *  @var _total_two_d_count: the total number of data points used to
         *        describe the field in the second dimension
         *  @var _total_three_d_count: the total number of data points used to
         *        describe the field in the third dimension
         */
	    const std::string _data_delim = "\t";
	    int _num_data_lines;
	    int _num_data_points_per_line;
    	int _total_one_d_count;
    	int _total_two_d_count;
    	int _total_three_d_count;
        ///@}

    	/** @name Pseudo Constructor Member Methods
         *  @brief Class wide member methods to be used within the class.
         *      These member methods are called and managed during the pseudo
         *      constructor processes.
	     */
        ///@{

        /** @name text_to_string_vector:
         *  @brief Converts a text file to a vector containing contents of each
         *      line and assigns this to the class wide variable _line_vector.
         *      Each element of the vector is a line of the text file in order.
         *      The Vector is of initially arbitrary size and scales to the
         *      size of the text file.
         *  @param field_file_name : the name of the text file containing the
         *      raw field information
         *  @var _line_vector: method wide container of the vector of strings.
         */
         void text_to_string_vector(std::string field_file_name);

         /** @name find_header_details:
          *  @brief Takes the header information (stored as strings in the @var
          *      _line_vector) and extracts the header details from them. These are stored as class wide variables (see @name Header Details).
          */
     	 void find_header_details();

         /** @name split_after_delim:
          *  @brief Takes a string and returns the substring after a provided
          *      delimiter.
          *  @param delim: a string denoting what is expected at the start of
          *      the main string.
          *  @param str_to_split: the main string to be split. This main string
          *      shall start with the delimiter and end with the important
          *      information.
          *  @return split_string: the string section containing the important
          *      information
          */
	     std::string split_after_delim(std::string delim, std::string str_to_split);

         /** @name string_to_bool:
          *  @brief Takes a string and returns a boolean. Can take capital or
          *      lowercase as well as 0 or 1.
          *  @param string: the string to be converted to a boolean
          *  @param delim: the message before the boolean that is being
          *      converted to help with error messages
          *  @return the corresponding boolean
          */
	     bool string_to_bool(std::string string, std::string delim);

         /** @name check_text_file_settings:
          *  @brief Checks the text file header details comply with what can be
          *      correctly interpreted here and sets up member variables (see
          *      @name Coordinate Sytem details). If there is no complicity,
          *      errors are thrown.
          *      Utilises @name is_expected_num_data_points_in_line submethod.
          */
         void check_text_file_settings();

         /** @name is_expected_num_data_points_in_line:
          *  @brief Checks that a line contains the expected number of data
          *      points. These data points must be a number or else they are
          *      not considered a data point.
          *  @param line: the line being checked for a certain number of data
          *      points.
          *  @param expected_num_data_points_per_line: the expected number of
          *      data points in the line.
          *  @return a boolean depending on both data points being numbers and
          *      the expected number of data points.
          */
         bool is_expected_num_data_points_in_line(std::string line, int expected_num_data_points_per_line);

         /** @name find_expected_number_data_points_per_line:
          *  @brief Finds the expected number of data points per line based on
          *      the stated file structure (how many dimensions are being
          *      described and whether the electric field is being defined or
          *       if the potential is being defined).
          *  @return the expected number of data points per data line.
          */
         int find_expected_number_data_points_per_line();

         /** @name get_values_from_data_string:
          *  @brief Extracts the values from a string containing data values
          *      and stores them as a vector of doubles where each element is a
          *      data value.
          *  @return line_values: the vector of data values stored in the line
          */
         std::vector<double> get_values_from_data_string(std::string line, int num_data_points_per_line);

         /** @name find_dust_potential_and_dust_radius:
          *  @brief Sets the dust potential @name _dust_potential as either
          *      that provided from DiMPl or from the first data point. If it
          *      is set from the first data point then a submethod is called to
          *      allow for conversion to necessary units.
          *      Furthermore, sets the dust radius @name _dust_radius method
          *      wide variable. Note for it to be set from the text file, then
          *      the text file itself cannot be normalised to Dust Radii.
          */
         void find_dust_potential_and_dust_radius(double dimpl_pars_dust_potential, double dimpl_pars_dust_radius);

         /** @name normalise_distance_to_dust_radius:
          *  @brief Normalises distances to the dust radius. All angles are
          *      expected to be in radians, but this could be changed here if
          *      necessary alongside an additional toggle in the header to
          *      allow for interpretation of other angular schemes (degrees
          *      etc.)
          *      Furthermore, should other distance normalisation schemes
          *      become useful in the future, then they should be normalised
          *      back to the dust radius at this juncture through the
          *      traditional use of a toggle in the text file header analagous
          *      to @name _is_pos_debye_length_normalised_truth.
          *  @param old_position_value: the unnormalised position value
          *      (typically from the text file).
          *  @param dimension: the dimension in which the position is valid
          *      (eg. the phi position in cylindrical coordinates has a
          *      dimension of 3).
          *  @return new_position_value: the normalised position value
          */
         double normalise_distance_to_dust_radius(double old_position_value, int dimension);

         /** @name normalise_potential_to_dimpl_scheme:
          *  @brief Normalises potentials to the scheme used by dimpl ( phi = (
          *      -eV)/(kBTe)).
          *      Currently only one potential normalisation scheme is
          *      permitted, but should alternative potential normalisation
          *      schemes be used in the future, then their translation to the
          *      standard dimpl normalisation scheme should be implemented here.
          *  @param old_potential: the unnormalised (or incorrectly normalised)
          *      potential value (typically from the text file).
          *  @return new_potential: the correctly normalised potential value
          */
         double normalise_potential_to_dimpl_scheme(double old_potential);

         /** @name normalise_E_Field_to_dimpl_scheme:
          *  @brief Normalises Electric fields to the scheme used by dimpl (see
          *      main.cpp, something like E_norm =(k_B T_e m/dust radius)*E_SI))
          *      Currently only one normalisation scheme is permitted, but
          *      should alternative normalisation schemes be used in the
          *      future, then their translation to the standard dimpl
          *      normalisation scheme should be implemented here.
          *  @param old_potential: the unnormalised (or incorrectly normalised)
          *      potential value (typically from the text file).
          *  @return new_potential: the correctly normalised potential value
          */
         double normalise_E_Field_to_dimpl_scheme(double old_EField);

         /** @name apply_theta_correction:
          *  @brief Modifies the theta value if necessary. This is since a theta
          *      value of exactly 0 or pi shall cause a badly defined phi in
          *      the spherical coordinate system. Modification performed by
          *      shifting thetas of 0 and pi by + and - delta respectively
          *  @param old_theta: the raw theta value provided (typically from the
          *      text file).
          *  @return new_theta: corrected theta value to be used
          */
         double apply_theta_correction(double old_theta);

         /** @name convert_data_line_to_multi_D_array:
          *  @brief Take raw data (in an ordered fashion) from the text file
          *      and transform it to be stored as a multi dimensional array.
          *  Procedure:
          *      Begins by recording the data dimensions in order to construct *      a 2d array.
          *      Fills 2d array with data from the text file. The data Points
          *      are normalised and corrected for dimpl interpretation.
          *      The 2d array is then converted to a multi dimensional array of
          *      field points.
          *      A count of how many points there are in each dimension is made.
          *      Then invalid Situations are accounted for (r<0; theta<0;
          *      theta>pi; phi<0; phi>2*pi; rho<0)
          *      Checks are made to ensure that there are enough field Points
          *      available if the user wishes to describe that dimension (min 3)
          *      Then a 3 dimensional grid (@var Field_Map_Matrix) is populated
          *      with the values, where the grid spans real space. Each point
          *      is a Field_Point object. Conversion from the 2D array to the
          *      3D grid is performed using the formula:
          *          (k*one_d_count*two_d_count) + j * one_d_count + i where k
          *      is the third dimension position, j is the second dimension
          *      position and i is the first dimension position. The
          *      Field_Point objects recognise a vector as being either
          *      cartesian, spherical or cylindrical depending on previously
          *      found type specified in file header.
          *  Note @var _two_d_array is stored as a class wide variable for
          *      memory issues with larger files on certain operating systems
          *     (stack vs heap) and emptied here since it is no longer required.
          */
         void convert_data_line_to_multi_D_array();

         /** @name partition_three_D_grid:
          *  @brief Used to handle the partitioning of the 3D grid of
          *      Field_Point objects. Works according to the number of
          *      dimensions being used to describe the field - it only
          *      partitions those dimensions that are described.
          *      If the dimension is not described, the `partition' is simply
          *      the single point.
          *      The partitions alluded to are the
          *      _(one/two/three)_d_pos_holder and
          *      _(one/two/three)_d_value_holder which are vector of vectors,
          *      in one direction the partition values, in the other direction
          *      additional granularity
          *  @submethod: find_per_dimension; actually does the partitioning
          */
	     void partition_three_D_grid();

         /** @name find_per_dimension:
          *  @brief Used to partition a Grid along a single dimension.
          *  Procedure:
          *      Specifies which dimension is being partitioned and manages
          *      that dimension of the 3D grid accordingly.
          *      Uses a `pos_holder' to account for the index of the grid point
          *      (the value 3 in a grid 2,3,4,5 shall have an index/pos of 1)
          *      and a `val_holder' to account for the corresponding value of
          *      interest (the position/distance in that dimension at the
          *      corresponding pos/index).
          *      The grid in that dimension is split up into sections (@param
          *      num_secs) where each section contains an additional level of
          *      granularity in that dimension.
          *      The first line is the two end points and the central points
          *      (0, n/2, n). For the pos, if n/2 is not an integer then it is
          *      rounded down. This is managed by the @submethod partition.
          *      The second line is the above line but with additional
          *      midpoints between all points (0, n/4, n/2, 3n/4, n). Again,
          *      non-integer pos's are rounded down.
          *      The third to penultimate (inclusive) line follow the same
          *      procedure (adding the midpoints each time)
          *      The final line is every single point.
          *      The corresponding class-wide member variable position and
          *      value holders for that dimension (
          *      _(one/two/three)_d_pos_holder and
          *      _(one/two/three)_d_value_holder) are updated accordingly
          *  @param dimension: int; the dimension in which the values are being
          *     partitioned
          */
     	 void find_per_dimension(int dimension);

         /** @name partition:
          *  @brief Used to find the midpoint between two positions.
          *  @param start: int; the index of the first position
          *  @param end: int; the index of the final position
          *  @var midway: int; the index of the midway position between the start and end. This will always return the rounded-down integer.
          *  @return vector containing the start, midaway and end index position
          */
	     std::vector<int> partition(int start, int end);

         /** @name find_val_from_field_point:
          *  @brief Used to get the corresponding value from a field point.
          *      The ordering of the cartesian, spherical and cylindrical are
          *      (x,y,z), (r, theta, phi) and (rho, z, phi) for dimension 1,2,3
          *      respectively.
          *  @param index: int; the index of the Field_Point being looked at
          *  @param dimension: int; the dimension being looked at
          *  @return value: double; the found corresponding value
          */
         double find_val_from_field_point(int index, int dimension);

         /** @name add_quadratic_fit_details_to_all_points:
          *  @breif Used to manage the finding of quadratic fits between points
          *      for all points in the Field_Map.
          *      Also calls a method on each Field_Point to record the fit
          *      details associated with the point.
          *      Furthermore calls a method on each Field_Point to find the
          *      electric field at that point.
          *      NOTE: the quadratic fits are only truly calculated for when
          *      the electric field is not defined, but this method should still
          *      be called when it is defined since the electric field at each
          *      field point must be found.
          *  @param pos: int[3]; the index of the position in the grid in terms
          *      of {first index position, second index position, third index
          *      position}
          *  @param fits: std::vector<std::vector<double>>; each line refers to
          *      one of the three dimensions and along each line is the a,b and
          *      c values referring to the fit v = a*s^2+b*s+c where v is the
          *      value of the field at that point (a potential) and s is the
          *      position in space along that dimension.
          *  @submethod find_fits
          */
 	     void add_quadratic_fit_details_to_all_points();

         /** @name find_fits:
          *  @brief Specifies the neighbouring points in each dimension
          *      required to find the fit.
          *      Procedure:
          *        Checks if the point being looked at is an edge point in
          *        which case it is treated differently
          *        Finds the neighbouring points (1 up and 1 down) in each of
          *        the three dimensions
          *        Retrieves the corresponding value at these points and the
          *        position in that dimension through the submethods @name get_specified_position and @name get_value
          *        Calls the @submethod calc_quad_fit to find the relationship
          *        of the central point to the neighbouring points for each
          *        dimension
          *    @return abc_s: the abc values for the quadratic fit for each of
          *        the three dimensions (stored as each row) with each of the
          *        three columns being a, b and c respectively
          */
     	 std::vector<std::vector<double>> find_fits(int position[3]);

         /** @brief calc_quad_fit:
          *    Calculates the quadratic fit equation between three points
          *    The fit is of the form y=ax^2+bx+c (quadratic)
          *    Procedure:
          *        Creates a set of three linear equations specifying the
          *        coordinates that a quadratic curve must fit through and
          *        solves them using a 3x3 matrix
          *        Checks if the three points are repeated etc. in which case
          *        the quadratic fit is not applicable and deals with these
          *        situations accordingly.
          *    @param p1: double[2]; the position (in the dimension of
          *        interest) and value (the potential) of the low point.
          *    @param p2: double[2]; the position (in the dimension of
          *        interest) and value (the potential) of the mid point.
          *    @param p3: double[2]; the position (in the dimension of
          *        interest) and value (the potential) of the high point.
          *    @return Y; the details of the fit written as {a,b,c}
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
         ///@}


         /** @name Main Member Methods.
          *  @brief Class wide member methods to be used within the class.
          *      These member methods are called and managed during the main
          *      loop of the programme (called many times) and so efficinecy is
          *      of upmost importance.
          */
         ///@{

         /** @name find_closest:
          *  @brief Uses the partitioned matrices of the index positions and
          *      values (the position in a certain dimension and coordinate
          *      system) to find the index of where the position being sought
          *      is.
          *    Procedure:
          *        Begin by checking if the position being sought is inside the
          *        Field map or not. If the position is outside of the field
          *        map then the edge position is said to be the closest and an
          *        instance of being outside the domain is recorded for later
          *        use.
          *        If it is inside the domain then begin by looking at the
          *        values in each line.
          *        For the first line, if the true value is above the midpoint
          *        then on the line beneath the next midpoint should be
          *        considered.
          *        This is enacted by keeping track of a @var `test_point'
          *        which is shifted accordingly for each subsequent line.
          *        This is continued up until the line before the penultimate
          *        line (since the final line is spaced differently)
          *        If the value is found to be below the `test_point' in the
          *        penultimate line, then the index below the `test_point'  is
          *        the base and the `test_point' is the upper. Conversely, if
          *        the value is found to be above the `test_point' in the
          *        penultimate line then the `test_point' is the base and the
          *        index above the `test_point'  is the upper.
          *        Next the spacing between the upper and base positions are
          *        found (@var num_to_check) and each of these position indices
          *        are checked on the final line to check which is the closest
          *        position and which is the second closest.
          *    @param value: the position being sought within the grid.
          *    @param pos_holder: the partitioned position holder matrix for
          *        the dimension in which the position value is being sought
          *    @param value_holder: the partitioned value holder matrix for the
          *        dimension in which the position value is being sought
          *    @int dimension: the dimension in which the position value is
          *        being sought
          *    @return closest_positions: the indices of the closest and second
          *        closest positions.
          */
	     std::vector<int> find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder, int dimension);
         ///@}


         /** @name Domain Checking Member Methods.
          *  @brief Private member methods to be used when checking the domain
          *      of the map.
          */
         ///@{
         /** @name check_coordinate_domain_symmetry.
          *  @brief Checks if the map domain is symmetrical along a certain
          *      dimension in a certain coordinate system and provides a
          *      warning message if it is not.
          *      Note: it is perhaps worth elevating this warning message to a
          *          formal error if this is found to be particularly bothersome
          *  @param dimension_to_check the dimension in which the symmetry is
          *      to be checked
          */
         void check_coordinate_domain_symmetry(int dimension_to_check);

         /** @name check_coordinate_domain_completeness.
          *  @brief Checks if the map domain is complete over the desired range
          *      Currently just checks that theta values, if provided, cover a
          *      range from 0 to pi and phi values, if provided, cover a range
          *      from 0 to 2 pi.
          *      Note: it is perhaps worth elevating this warning message to a
          *          formal error if this is found to be particularly bothersome
          *  @param dimension_to_check the dimension in which the symmetry is
          *      to be checked
          */
         void check_coordinate_domain_completeness(int dimension_to_check);

         /** @name find_domain_limits.
          *  @brief Finds the critical limits of the map domain and records
          *      them to the class wide variables @var _map_z_min, @var
          *      _map_z_max @var _map_rho_max
          *      For the situation where the map is described in spherical
          *      coordinates, but cylindrical injection is used, the largest
          *      cylinder to fit inside the sphere is used.
          */
         void find_domain_limits();

         /** @name check_dust_grain_map_coverage.
          *  @brief Checks that the map covers a domain that is sufficiently
          *      close to the dust grain surface (otherwise large errors are
          *       possible)
          *      This is only performed for the spherical and cylindrical case
          *      since the Cartesian case should already be large enough and
          *      symmetrical (these checks are performed elsewhere)
          */
         void check_dust_grain_map_coverage();

         /** @name check_domain.
          *  @brief Manages the calling of the various submethods that check
          *      the domain described by the map is acceptable. Currently only
          *      warning messages are provided for unacceptable domains but it
          *      is perhaps worth elevating this warning message to a formal
          *      error if this is found to be particularly bothersome.
          *  Checks performed:
          *     Check for symmetry and completeness.
          *     Check for logical choice of coordinate system.
          *     Check the z and radial injection limits set by dimpl are
          *         covered by the custom map
          *     Check for dust limits
          */
         void check_domain();

         /** @name get_shrunk_point
          *  @brief Finds the nearby point in the map that will encapsulate the
          *      new value. This is to be used when seeking to shrink the map.
          *  @param new_val: the new position value at the limit of the map
          *  @param old_val: the old limit of the map value
          *  @param dimension: the dimension in which the position values exist
          *  @param are_values_positive: whether the new map limiting value is
          *      positive or negative
          *  @return new_pos: the position in the old map in which the new map
          *      limit exists
          */
         int get_shrunk_point(double new_val, double old_val, int dimension, bool are_values_positive);
         ///@}

         /** @name General Member Methods.
          *  @brief Private member methods to be used generally throughout the
          *      class
          */
         ///@{

         /** @name get_coordinate_name.
          *  @brief Gets the name that a certain dimension represents in the
          *      used coordinate system.
          *  @param dimension: the dimension for which the name is required
          *  @return coordinate_name: the name (as a string) of this
          *      dimension in this coordinate system
          */
         std::string get_coordinate_name(int dimension);

         /** @name check_if_vectors_differ.
          *  @brief Goes through elements of two vectors to see if they differ
          *  @param one: the first vector to check
          *  @param two: the second vector to check
          *  @return is_different: whether the vectors differ or not
          */
         bool check_if_vectors_differ(std::vector<int> one, std::vector<int> two);
         ///@}

    public:
	/** @name Constructors
	 *  @brief constructs a Field_Map class
	 */
	///@{
    /** @name Default Constructor.
     *  @brief Default constructor. Creates a field map object shell when called
     *      but this is unuseable until the pseudo constructor
     *      (construct_in_full method) is subsequently called to enable full
     *      construction.
     */
	Field_Map();
	/** @name construct_in_full
     *  @brief Pseudo Constructor. Manages entire setup of a Field_Map object
     *      through the calling of various methods each of which is essential
     *      for the full description of a Field_Map.
     *      Procedure:
     *          Call method to transform text file to vector of strings
     *          Call method to find the header details for the Field Map
     *              construction
     *          Call method to check the details provided in the text file
     *              header are acceptable
     *          Call method to assign to find the Debye Length and Dust Radius
     *          Call method to transform vector of strings to multi dimensional
     *              array of Field Point objects
     *          Call method to partition the multi dimensional array for easy
     *              position searching
     *          Calls method to find the relationship between neighbouring
     *              Field_Points to create a quadratic fit between them.
     *          Call method to check that the extent of the map domain complies
     *              with what this software can understand
     *          Call method to find the extent of the map domain for subsequent
     *              checking with DiMPl pars.
     *   @param field_file_name: the name of the file containing the field
     *          details. This should be stored inside the Custom Fields
     *          directory.
     *  @param dimpl_pars_dust_potential: the dust potential provided by the
     *          DiMPl pars input.
     *  @param dimpl_pars_dust_radius: the dust radius provided by the dimpl
     *         pars input.
     *  @param dimpl_pars_debye_length: the Debye Length as calculated from the
     *         DiMPl pars input as found in the main.cpp file.
     *  @param dimpl_pars_eTemp: the electron temperature provided by the dimpl
     *         pars input.
     *  @param is_spherical_injection: whether spherical injection is to be
     *         used by DimPl or not.
     *  @param is_cylindrical_injection: whether cylindrical injection is to be
     *         used by DimPl or not.
     *  @param a1: the dust radius 1/(semix^2) paramater used by DiMPl
     *  @param a2: the dust radius 1/(semiy^2) paramater used by DiMPl
     *  @param a3: the dust radius 1/(semiz^2) paramater used by DiMPl
     */
	void construct_in_full(std::string field_file_name, double dimpl_pars_dust_potential, double dimpl_pars_dust_radius, double dimpl_pars_debye_length, double dimpl_pars_eTemp, bool is_spherical_injection, bool is_cylindrical_injection, double a1, double a2, double a3);
    ///@}

    /** @name Public Main Member Methods.
     *  @brief Class wide member methods to be used within the class.
     *      These member methods are called and managed during the main
     *      loop of the programme (called many times) and so efficinecy is
     *      of upmost importance.
     */
    ///@{
    /** @name find_approx_value:
     *  @brief Searches in the Field_Map for the approximate field value at the
     *      position provided. The position provided does not need to be an
     *      exact point described in the Field_Map.
     *      The @submethod find_closest is used to find all of the nearby field
     *      points (number required depends on how many dimensions the map is
     *      being described in).
     *      The electric field is found according to each field point (if the
     *      field is to be found on the curve then the abc details are used and
     *      it is found in the Field_Point class) otherwise the electric field
     *      at that point is used.
     *      The distance between the nearby field points and desired point is
     *      found and the electric field according to each point is weighted
     *      according to the inverse square of this distance
     *      Checks are performed on how to manage matters should the desired
     *      point lie exactly on a field point
     *  @param point: the position threevector of the point at which the
     *      electric field value is required.
     *  @return approximate_E_Field_val; approximate electric field value at
     *      the point provided
     */
	threevector find_approx_value(threevector point);
    ///@}

    /** @name Public Domain Checking Member Methods.
     *  @brief Class wide member methods to be used within the class.
     *      These member methods are called when checking the domain limits of
     *      the field map.
     */
    ///@{
    /** @name set_new_map_limits.
     *  @brief Sets the critical limits of the map domain when a rho limit is
     *      provided (from the cylindrical injection) in order to maximise the
     *      size of the cylinder that would fit inside the sphere.
     *  @param new_rho_limit the radius of the cylinder being used for
     *      cylindrical injection
     */
    void set_new_map_limits(double new_rho_limit);

    /** @name shrink_map_to_fit.
     *  @brief Shrinks the map to fit the domain size specified by the
     *      injection in dimpl.
     *      Uses the @submethod get_shrunk_point to find the points in the old
     *      map that would encapsulate the new limits.
     *      If the new position limits differ from the old position limits then
     *      shrinking is required
     *      Shrinking is preformed by spanning the old map between the new
     *      limits in each dimension and recording it to a new map. This map is
     *      then told to replace the old map and the various class members are
     *      updated accordingly
     *  @param z_min the minimum z value of the new domain
     *  @param z_max the maximum z value of the new domain
     *  @param rho_max the maximum rho value of the new domain
     */
    void shrink_map_to_fit(double z_min, double z_max, double rho);

    /** @name check_num_times_outside_domain.
     *  @brief Checks the number of times there were instances of seeking values
     *     outside of the domain and prints the number of times this occurred if
     *     it is non-zero
     */
    void check_num_times_outside_domain();
    ///@}

    /** @name get methods
     *  @brief returns the various respective class members
     */
    ///@{
	inline double get_map_z_max(){return _map_z_max;};
	inline double get_map_z_min(){return _map_z_min;};
	inline double get_map_rho_max(){return _map_rho_max;};
	inline double get_debye_length(){return _debye_length;};
	inline double get_dust_radius(){return _dust_radius;};
	inline double get_dust_potential(){return _dust_potential;};
    ///@}

}; //end of class

#endif /*__FIELD_MAP_H_INCLUDED__ */
