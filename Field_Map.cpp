/** @file Field_Map.cpp
 *  @brief Implementation of Field_Map class for use in describing Physics fields
 *
 *  Field_Map Objects contain details of a field comprised of many Field_Point objects
 *
 *  @author Daniel Greenhouse (dg2020@ic.ac.uk (until 01.09.2020), dmgreenhouse@outlook.com (thereafter))
 *      Written immitating the style of main DiMPl author Luke Simons (ls5115@ic.ac.uk)
 *
 *  @bugs No known bugs.
 */

#include "Field_Map.h" //!< For the Field_Map class hereby implemented

/** Constructors */
Field_Map::Field_Map(){}

void Field_Map::construct_in_full(std::string field_file_name, double dimpl_pars_dust_potential, double dimpl_pars_dust_radius, double dimpl_pars_debye_length, double dimpl_pars_eTemp, bool is_spherical_injection, bool is_cylindrical_injection, double a1, double a2, double a3)
{
    // assign DiMPl pars to this class
    _eTemp = dimpl_pars_eTemp;
    _debye_length = dimpl_pars_debye_length;
    _is_spherical_injection = is_spherical_injection;
    _is_cylindrical_injection = is_cylindrical_injection;
    _a1 = a1;
    _a2 = a2;
    _a3 = a3;
    // make subsequent calls to methods for full map construction
    if(_DEBUG){ std::cout<<"Opened constructor and assigned Parameters."<<std::endl;}

    text_to_string_vector(field_file_name);
    if(_DEBUG){std::cout<<"text to string vector complete."<<std::endl;}

    find_header_details();
    if(_DEBUG){std::cout<<"Header details found."<<std::endl;}

    check_text_file_settings();
    if(_DEBUG){std::cout<<"Checked Text File settings."<<std::endl;}

    find_dust_potential_and_dust_radius(dimpl_pars_dust_potential, dimpl_pars_dust_radius);
    if(_DEBUG){std::cout<<"Found the Debye Length and Dust Radius."<<std::endl;}

    convert_data_line_to_multi_D_array();
    if(_DEBUG){std::cout<<"Convereted data to multi d array."<<std::endl;}

    partition_three_D_grid();
    if(_DEBUG){std::cout<<"Partitioned three D grid."<<std::endl;}

    add_quadratic_fit_details_to_all_points();
    if(_DEBUG){std::cout<<"Fitted all points with a quadratic fit."<<std::endl;}

    check_domain();
    if(_DEBUG){std::cout<<"Checked domain."<< std::endl;}

    find_domain_limits();
    if(_DEBUG){std::cout<<"Found domain limits."<<std::endl;}
}


void Field_Map::text_to_string_vector(std::string field_file_name)
{
    std::ifstream in(field_file_name);
    std::string single_line;
    if (in){
        while(getline(in, single_line)) {
	        _line_vector.push_back(single_line);
        }
    }
}

void Field_Map::find_header_details()
{
    _summary_line = split_after_delim(_summary_delim, _line_vector[0]);
    _dimension_val = std::stoi(split_after_delim(_dimension_delim, _line_vector[1])); // and convert string to integer
    _coord_type = split_after_delim(_coord_delim, _line_vector[2]);
    _ordered_truth = string_to_bool(split_after_delim(_ordered_delim, _line_vector[3]), _ordered_delim); // and convert string to boolean
    _is_E_Field_defined_truth = string_to_bool(split_after_delim(_is_E_Field_defined_delim, _line_vector[4]), _is_E_Field_defined_delim); // and convert string to boolean
    _is_pos_dust_radius_normalised_truth = string_to_bool(split_after_delim(_is_pos_dust_radius_normalised_delim, _line_vector[5]), _is_pos_dust_radius_normalised_delim); // and convert string to boolean
    _is_pos_debye_length_normalised_truth = string_to_bool(split_after_delim(_is_pos_debye_length_normalised_delim, _line_vector[6]), _is_pos_debye_length_normalised_delim); // and convert string to boolean
    _is_potential_normalised_truth = string_to_bool(split_after_delim(_is_potential_normalised_delim, _line_vector[7]), _is_potential_normalised_delim); // and convert string to boolean
    _is_initial_dust_pot_truth = string_to_bool(split_after_delim(_is_initial_dust_pot_delim, _line_vector[8]), _is_initial_dust_pot_delim); // and convert string to boolean
    _is_initial_dust_rad_truth = string_to_bool(split_after_delim(_is_initial_dust_rad_delim, _line_vector[9]), _is_initial_dust_rad_delim); // and convert string to boolean
}


std::string Field_Map::split_after_delim(std::string delim, std::string str_to_split)
{
    const int end_delim_pos = str_to_split.find(_end_delim);
    const int delim_length = delim.length();
    const int delim_pos = str_to_split.find(delim);
    const int start_index = delim_pos + delim_length;
    const int str_to_split_length = str_to_split.length();
    if (end_delim_pos<=0 || delim_pos<0 || start_index>=str_to_split_length){
	std::cerr<<"ERROR: the header line has not been correctly interpreted."
	<<"\n\tCheck that the header delimiters are correctly spelled."
	<<"\n\t\tExpected header delimiter: '"<<delim<<"'"
	<<"\n\t\tTrying to extract information from line: '"<<str_to_split<<"'"
	<<"\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	throw std::exception();
    }
    const std::string split_string = str_to_split.substr(start_index, end_delim_pos-start_index);
    return split_string;
}

bool Field_Map::string_to_bool(std::string string, std::string delim){
    const std::string unknown_boolean_warning = "WARNING: the boolean has not been recognised.\n\tPlease correct the Custom Field Text File and see the README for help.";
    bool string_as_boolean = true; // initially assume ordered
    if (string.size() == 1){
        // the number case
	    if (string=="0"){
	        string_as_boolean = false;
	    }
	    else if (string=="1"){
	        string_as_boolean = true;
	    }
	    else {
	        std::cerr<< unknown_boolean_warning <<std::endl;
	        string_as_boolean = true;
		    throw std::exception();
	    }
    } else {
        // the text boolean case
        std::transform(string.begin(), string.end(), string.begin(), [](unsigned char c) {return std::tolower(c); }); // convert string to lower case
	    if (string=="true"){string_as_boolean = true;}
	    else if (string=="false"){string_as_boolean = false;}
        else {
	        std::cerr << unknown_boolean_warning << "\tDelimiter stated to be: "<<delim<<"\tString stated to be: "<<string<<std::endl;
	        string_as_boolean = true;
		    throw std::exception();
	    }
    }
    return string_as_boolean;
}

void Field_Map::check_text_file_settings(){
    // Check the overriding coordinate system.
    const std::string dimension_mismatch_warning = "ERROR: coordinate system dimensions do not match those specified under `Number of Dimensions'.\n\tPlease correct the Custom Field Text File and see the README for help.";
    if (_coord_type == _CART_1 || _coord_type == _CART_2 || _coord_type == _CART_3){
        _is_cartesian = true;
        if (_coord_type == _CART_1 && _dimension_val!=1){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        } else if (_coord_type == _CART_2 && _dimension_val!=2){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        } else if (_coord_type == _CART_3 && _dimension_val!=3){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        }
    } else if (_coord_type == _SPHER_1 || _coord_type == _SPHER_2 || _coord_type == _SPHER_3){
        _is_spherical = true;
        if (_coord_type == _SPHER_1 && _dimension_val!=1){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        } else if (_coord_type == _SPHER_2 && _dimension_val!=2){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        } else if (_coord_type == _SPHER_3 && _dimension_val!=3){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        }
    } else if (_coord_type == _CYL_1 || _coord_type == _CYL_2 || _coord_type == _CYL_3){
        _is_cylindrical = true;
        if (_coord_type == _CYL_1 && _dimension_val!=1){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        } else if (_coord_type == _CYL_2 && _dimension_val!=2){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        } else if (_coord_type == _CYL_3 && _dimension_val!=3){
            std::cerr<<dimension_mismatch_warning<<std::endl;
	    throw std::exception();
        }
    } else {
        std::cerr<<"ERROR: coordinate system is not recognised.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }
    if (!_ordered_truth){
	    std::cerr<<"ERROR: method only currently works for ordered files.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }
    // Check that the data is separated by tabs
    const std::string expected_data_line = _line_vector[_num_lines_in_header];
    const int expected_num_data_points_per_line = find_expected_number_data_points_per_line();
    const bool is_expected_data_line_ok = is_expected_num_data_points_in_line(expected_data_line, expected_num_data_points_per_line);
    if (!is_expected_data_line_ok){
	    std::cerr<<"ERROR: the data in the custom file has not been correctly delimited or has an additional line in the header. Tabs should be used between values.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl << expected_data_line << std::endl;
	    throw std::exception();
    }
    // Check that the data starts where it is expected to start
    const std::string unexpected_data_line = _line_vector[_num_lines_in_header-1];
    const bool is_unexpected_data_line_ok = is_expected_num_data_points_in_line(unexpected_data_line, expected_num_data_points_per_line);
    if (is_unexpected_data_line_ok){
	    std::cerr<<"ERROR: the data in the custom file contains an in-complete header. \n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }

    // Check if there is some valid normalisation for the position values
    if (!_is_pos_dust_radius_normalised_truth && !_is_pos_debye_length_normalised_truth){
        std::cerr<<"ERROR: a valid normalisation scheme for the position data is required. \n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }
    else if (_is_pos_dust_radius_normalised_truth && _is_pos_debye_length_normalised_truth){
        std::cerr<<"ERROR: the position data must be normalised according to one scheme only. \n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }
    // Check if there is some valid normalisation for the potential values
    if (!_is_potential_normalised_truth){
	    if (_is_E_Field_defined_truth){
	        std::cerr<<"ERROR: a valid normalisation scheme for the electric field values is required. \n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
	    } else {
	        std::cerr<<"ERROR: a valid normalisation scheme for the potential values is required. \n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
	    }
    }
    // Check if the potential is allowed to be found from the text file
    if (_is_initial_dust_pot_truth){
	    if (_is_E_Field_defined_truth){
	        std::cerr<<"ERROR: the dust grain potential cannot be extracted from the Custom Field Map File when the electric field is being defined.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
	    }
	    if (!_is_spherical) {
	        std::cerr<<"ERROR: the dust grain potential cannot be extracted from the Custom Field Map File when the geometry is not spherical.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
	    }
    }
    // Check if the dust radius is allowed to be found from the text file
    if (_is_initial_dust_rad_truth){
	    if (!_is_spherical) {
	        std::cerr<<"ERROR: the dust grain radius cannot be extracted from the Custom Field Map File when the geometry is not spherical.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
	    }
	    if (_is_pos_dust_radius_normalised_truth) {
	        std::cerr<<"ERROR: the dust grain radius cannot be extracted from the Custom Field Map File since the Field Map values are already normalised to the dust grain radius.\n\tTo proceed, enter the dust grain radius to DiMPl and set `Find Dust Grain Potential from this custom field map' to false in the custom field map text file.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
	    }
    }
}

bool Field_Map::is_expected_num_data_points_in_line(std::string line, int expected_num_data_points_per_line){
    bool is_expected_num_data_points= true;
    int true_data_points_count = 0;
    bool is_value_a_number = true;
    bool is_at_final = false;
    while (!is_at_final){
        int pos = line.find(_data_delim);
	int comparison_string_pos = std::string::npos;
	    if (pos != comparison_string_pos){
	        std::string line_data = line.substr(0, pos);
	        if (!line_data.empty() && line_data.find_first_not_of("0123456789-+eE.")!=std::string::npos){
		        is_value_a_number = false;
	        } else {
		        true_data_points_count += 1;
	        }
	        line.erase(0, pos + _data_delim.length());
	    } else {
	        if (!line.empty()&&line.find_first_not_of("0123456789-+eE.")!=0){
		        true_data_points_count += 1;
	        }
	        is_at_final = true;
	    }
    }
    if (true_data_points_count != expected_num_data_points_per_line){
        is_expected_num_data_points = false;
    }
    return is_expected_num_data_points && is_value_a_number;
}

int Field_Map::find_expected_number_data_points_per_line(){
    int num_data_points_per_line;
    if (_is_E_Field_defined_truth) {
	    num_data_points_per_line=2*_dimension_val;
    }
    else {
	    num_data_points_per_line=1+_dimension_val;
    }
    return num_data_points_per_line;
}

std::vector<double> Field_Map::get_values_from_data_string(std::string line, int num_data_points_per_line){
    std::string line_data[num_data_points_per_line];
    std::vector<double> line_values(num_data_points_per_line);
    for (int j=0; j<num_data_points_per_line; j++){
	    int pos = line.find(_data_delim);
	    line_data[j] = line.substr(0, pos);
	    line.erase(0, pos + _data_delim.length());
        std::stringstream value(line_data[j]);

	    if (value.fail()) {
            std::string error_string = "Unable to format ";
            error_string += line_data[j].c_str();
            error_string += " as a number.";
            std::cerr<<"ERROR: a value in the Custom Field Map file cannot be interpreted."
	        <<"\n\t"+error_string
	        <<"\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	        throw std::exception();
        } else {
            value >> line_values[j];
	    }
    }
    return line_values;
}

void Field_Map::find_dust_potential_and_dust_radius(double dimpl_pars_dust_potential, double dimpl_pars_dust_radius){
    // get the first data line from the file
    std::string first_line = _line_vector[_num_lines_in_header];
    const int num_data_points_per_line = find_expected_number_data_points_per_line();
    std::vector<double> first_line_values = get_values_from_data_string(first_line, num_data_points_per_line);

    if(_is_initial_dust_pot_truth){
	    _dust_potential = normalise_potential_to_dimpl_scheme(first_line_values[0]);
    } else {
	    _dust_potential = dimpl_pars_dust_potential;
    }
    if(_is_initial_dust_rad_truth){
	    double text_file_normalised_dust_radius;
	    if(_is_E_Field_defined_truth){
	        text_file_normalised_dust_radius = first_line_values[1];
	    } else{
	        text_file_normalised_dust_radius = first_line_values[_dimension_val];
	    }
	    if (_is_pos_debye_length_normalised_truth){
	        _dust_radius = text_file_normalised_dust_radius*_debye_length;
        } else if (_is_pos_dust_radius_normalised_truth) {
	        std::cerr<<"WARNING! The text file is normalised to dust radius so the dust radius cannot be correctly determined from it."<<std::endl;
            _dust_radius = dimpl_pars_dust_radius;
            throw std::exception();
        }
    } else {
	    _dust_radius = dimpl_pars_dust_radius;
    }
}

double Field_Map::normalise_distance_to_dust_radius(double old_position_value, int dimension){
    // If angular, assume in radians
    double new_position_value;
    if ( (_is_spherical && dimension == 2) || (_is_spherical && dimension == 3 ) || (_is_cylindrical && dimension == 3)){
        new_position_value = old_position_value;
    }
    // If not angular, normalise if not already in normalised scheme
    else {
        if (_is_pos_debye_length_normalised_truth){
	        new_position_value = old_position_value*_debye_length/_dust_radius;
        }
        else if (_is_pos_dust_radius_normalised_truth) {
            new_position_value = old_position_value;
        }
    }
    return new_position_value;
}

double Field_Map::normalise_potential_to_dimpl_scheme(double old_potential){
    double new_potential;
    if (_is_potential_normalised_truth){
	    new_potential = -old_potential;
    }
    // Currently this is the only normalisation scheme allowed, so below should not be called
    else {
        std::cerr<<"WARNING! The text file uses an unknown potential normalisation scheme. Please implement the correct method in normalise_potential_to_dimpl_scheme of Field_Map.cpp. Proceeding using standard normalisation scheme for now."<<std::endl;
        throw std::exception();
        new_potential = old_potential;
    }
    return new_potential;
}

double Field_Map::normalise_E_Field_to_dimpl_scheme(double old_EField){
    double new_EField;
    if (_is_potential_normalised_truth){
	    new_EField = -old_EField;
    }
    // Currently this is the only normalisation scheme allowed, so below should not be called
    else {
        std::cerr<<"WARNING! The text file uses an unknown electric field normalisation scheme. Please implement the correct method in normalise_E_Field_to_dimpl_scheme of Field_Map.cpp. Proceeding using standard normalisation scheme for now."<<std::endl;
        throw std::exception();
        new_EField = old_EField;
    }
    return new_EField;
}

double Field_Map::apply_theta_correction(double old_theta){
    double new_theta = old_theta;
    const double delta = 1e-15;
    if (_is_spherical) {
	    if (old_theta==0.0){
            new_theta = delta;
	    } else if (old_theta==PI) {
            new_theta -= delta;
        }
    }
    return new_theta;
}

void Field_Map::convert_data_line_to_multi_D_array(){
    // Begin by finding the size
    double two_d_array_width;
    const int num_data_points_per_line = find_expected_number_data_points_per_line();
    if (_is_E_Field_defined_truth) {two_d_array_width=6;}
    else {two_d_array_width=4;}
    _num_data_lines = _line_vector.size() - _num_lines_in_header;
    // Create array of correct shape
    _two_D_array = std::vector<std::vector<double>>(two_d_array_width, std::vector<double>(_num_data_lines));
    // Fill array with data points from text file.
    // Ensure data points are normalised to necessary scheme
    for (int i=0; i<_num_data_lines; i++){
	    std::string line = _line_vector[i+_num_lines_in_header];
	    std::vector<double> line_values = get_values_from_data_string(line, num_data_points_per_line);
        if (_dimension_val==1){
            if (_is_E_Field_defined_truth){
		        _two_D_array[0][i] = normalise_E_Field_to_dimpl_scheme(line_values[0]);
                _two_D_array[1][i] = normalise_E_Field_to_dimpl_scheme(apply_theta_correction(0.0));
                _two_D_array[2][i] = normalise_E_Field_to_dimpl_scheme(0.0);
		        _two_D_array[3][i] = normalise_distance_to_dust_radius(line_values[1],1);
		        _two_D_array[4][i] = normalise_distance_to_dust_radius(apply_theta_correction(0.0),2);
                _two_D_array[5][i] = normalise_distance_to_dust_radius(0.0,3);
            }
            else {
		        _two_D_array[0][i] = normalise_potential_to_dimpl_scheme(line_values[0]);
		        _two_D_array[1][i] = normalise_distance_to_dust_radius(line_values[1],1);
		        _two_D_array[2][i] = normalise_distance_to_dust_radius(apply_theta_correction(0.0),2);
                _two_D_array[3][i] = normalise_distance_to_dust_radius(0.0,3);
            }
        }
        if (_dimension_val==2){
            if (_is_E_Field_defined_truth){
		        _two_D_array[0][i] = normalise_E_Field_to_dimpl_scheme(line_values[0]);
		        _two_D_array[1][i] = normalise_E_Field_to_dimpl_scheme(apply_theta_correction(line_values[1]));
                _two_D_array[2][i] = normalise_E_Field_to_dimpl_scheme(0.0);
		        _two_D_array[3][i] = normalise_distance_to_dust_radius(line_values[2], 1);
		        _two_D_array[4][i] = normalise_distance_to_dust_radius(apply_theta_correction(line_values[3]),2);
                _two_D_array[5][i] = normalise_distance_to_dust_radius(0.0,3);
            }
            else {
		        _two_D_array[0][i] = normalise_potential_to_dimpl_scheme(line_values[0]);
		        _two_D_array[1][i] = normalise_distance_to_dust_radius(line_values[1],1);
		        _two_D_array[2][i] = normalise_distance_to_dust_radius(apply_theta_correction(line_values[2]),2);
                _two_D_array[3][i] = normalise_distance_to_dust_radius(0.0,3);
            }
        }
        if (_dimension_val==3){
            if (_is_E_Field_defined_truth){
		        _two_D_array[0][i] = normalise_E_Field_to_dimpl_scheme(line_values[0]);
		        _two_D_array[1][i] = normalise_E_Field_to_dimpl_scheme(apply_theta_correction(line_values[1]));
                _two_D_array[2][i] = normalise_E_Field_to_dimpl_scheme(line_values[2]);
		        _two_D_array[3][i] = normalise_distance_to_dust_radius(line_values[3],1);
		        _two_D_array[4][i] = normalise_distance_to_dust_radius(apply_theta_correction(line_values[4]),2);
                _two_D_array[5][i] = normalise_distance_to_dust_radius(line_values[5],3);
            }
            else {
		        _two_D_array[0][i] = normalise_potential_to_dimpl_scheme(line_values[0]);
		        _two_D_array[1][i] = normalise_distance_to_dust_radius(line_values[1],1);
		        _two_D_array[2][i] = normalise_distance_to_dust_radius(apply_theta_correction(line_values[2]),2);
                _two_D_array[3][i] = normalise_distance_to_dust_radius(line_values[3],3);
            }
        }
    }


    // Convert grid representation to multi dimensional array of Field Points
    // Begin by finding number of points in each dimension
    int one_d_count = 1;
    int two_d_count = 1;
    int three_d_count = 1;
    int first_dimension_pos = 1;
    if (_is_E_Field_defined_truth){first_dimension_pos = 3;}
    double memory_holder = _two_D_array[first_dimension_pos][0]; // start off with the first value in the first dimension

    for (int i=1; i< _num_data_lines+1; i++){
        const double first_d_current_val = _two_D_array[first_dimension_pos][i];
        if ( first_d_current_val < memory_holder){
	        break;
	    }
	    else {
	        memory_holder = first_d_current_val;
	        one_d_count+=1;
	    }
    }
    if (_dimension_val==2||_dimension_val==3){
        memory_holder = _two_D_array[first_dimension_pos+1][0]; // start off with the first value
        for (int i=1; i<_num_data_lines/one_d_count+1; i++){
            const int step = i*one_d_count;
            const double second_d_current_val = _two_D_array[first_dimension_pos+1][step];
            if (second_d_current_val<memory_holder){break;}
	        else {
	            memory_holder = second_d_current_val;
	            two_d_count +=1;
	        }
        }
        if (_dimension_val==3){
        memory_holder = _two_D_array[first_dimension_pos+2][0]; // start off with the first value
            for (int i=1; i<_num_data_lines/(one_d_count*two_d_count)+1; i++){
                const int step = i*one_d_count*two_d_count;
		if (step > _num_data_lines){break;}
                const double third_d_current_val = _two_D_array[first_dimension_pos+2][step];
                if (third_d_current_val<memory_holder){
                    break;
                }
    	        else {
    	            memory_holder = third_d_current_val;
    	            three_d_count +=1;
    	        }
            }
        }
    }

    // Sort Out Invalid Situations:

    //    1. r is less than or equal to 0
    //    2. theta is less than 0 or greater than pi
    //    3. phi is less than 0 or greater than two pi
    //    4. rho is less than or equal to 0

    int one_d_low_error_count = 0; // initialise by saying there is no error
    int two_d_low_error_count = 0; // initialise by saying there is no error
    int two_d_high_error_count = 0; // initialise by saying there is no error
    int three_d_low_error_count = 0; // initialise by saying there is no error
    int three_d_high_error_count = 0; // initialise by saying there is no error
    if (_is_spherical || _is_cylindrical){
        // Make a Memory of position values in each dimension
        for (int i=0; i<one_d_count; i++){
	        if (_two_D_array[first_dimension_pos][i] <= 0.0){
	            one_d_low_error_count+=1;
	        }
	    }
       	if (_dimension_val ==2||_dimension_val==3){
	        const double TWO_PI = 2*PI;
            if (_is_spherical){
	            for (int i=0; i<two_d_count; i++){
	                double two_d_pos = _two_D_array[first_dimension_pos+1][i*one_d_count];
		            if (two_d_pos <= 0.0){
		                two_d_low_error_count+=1;
		            } else if (two_d_pos >=PI){
		                two_d_high_error_count += 1;
		            }
	            }
	        }
	        if (_dimension_val==3){
                for (int i=0; i<three_d_count; i++){
	                double three_d_pos = _two_D_array[first_dimension_pos+2][i*one_d_count*two_d_count];
		            if (three_d_pos < 0.0){
		                three_d_low_error_count+=1;
		            } else if (three_d_pos >TWO_PI){
		                three_d_high_error_count += 1;
		            }
	            }
	        }
        }
	    if (one_d_low_error_count != 0){
	        std::cout<<"WARNING: There was a field defined at r/rho values below 0. These have been removed but caution should be taken."<<std::endl;
	    }
	    if (two_d_low_error_count != 0){
	        std::cout<<"WARNING: There was a field defined at theta values below 0. These have been removed but caution should be taken."<<std::endl;
	    }
	    if (two_d_high_error_count != 0){
	        std::cout<<"WARNING: There was a field defined at theta values above pi. These have been removed but caution should be taken."<<std::endl;
	    }
	    if (three_d_low_error_count != 0){
	        std::cout<<"WARNING: There was a field defined at phi values below 0. These have been removed but caution should be taken."<<std::endl;
	    }
	    if (three_d_high_error_count != 0){
	        std::cout<<"WARNING: There was a field defined at phi values above 2pi. These have been removed but caution should be taken."<<std::endl;
	    }
    }

    _total_one_d_count = one_d_count-one_d_low_error_count;
    _total_two_d_count = two_d_count-two_d_low_error_count-two_d_high_error_count;
    _total_three_d_count = three_d_count - three_d_low_error_count - three_d_high_error_count;

    if (_total_one_d_count<3){
	    std::cerr<<"ERROR: the custom field does not contain enough valid values in the first dimension.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }

    if (((_dimension_val==2||_dimension_val==3) && _total_two_d_count<3)||(_dimension_val==1 && _total_two_d_count<=0)){
	    std::cerr<<"ERROR: the custom field does not contain enough valid values in the second dimension.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }

    if ((_dimension_val==3 && _total_three_d_count<3)||((_dimension_val==1||_dimension_val==2) && _total_three_d_count<=0)){
	    std::cerr<<"ERROR: the custom field does not contain enough valid values in the third dimension.\n\tPlease correct the Custom Field Text File and see the README for help."<<std::endl;
	    throw std::exception();
    }

    // Populate multi dimensional array
    std::vector<std::vector<std::vector<Field_Point>>> Field_Map_Matrix(one_d_count, std::vector<std::vector<Field_Point>>(two_d_count, std::vector<Field_Point>(three_d_count)));
    int i_count = 0;
    int j_count = 0;
    int k_count = 0;
    for (int i_true=one_d_low_error_count; i_true<one_d_count; i_true++){
        for (int j_true=two_d_low_error_count; j_true<two_d_count-two_d_high_error_count; j_true++) {
	        for (int k_true=three_d_low_error_count; k_true<three_d_count-three_d_high_error_count; k_true++){
	            int two_d_array_pos = (k_true*one_d_count*two_d_count) + j_true * one_d_count + i_true; // find corresponding position in the two d grid
                if (_is_cartesian){
                    if (_is_E_Field_defined_truth) {
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], _two_D_array[1][two_d_array_pos], _two_D_array[2][two_d_array_pos], threevector(_two_D_array[3][two_d_array_pos], _two_D_array[4][two_d_array_pos], _two_D_array[5][two_d_array_pos]), _is_cartesian, _is_spherical, _is_cylindrical);
                    } else {
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], threevector(_two_D_array[1][two_d_array_pos], _two_D_array[2][two_d_array_pos], _two_D_array[3][two_d_array_pos]));
                    }
                }
                else if (_is_spherical) {
		        // Note spherical defined in text file as r,theta,phi
		        // Threevector in order: r, theta, phi
                    if (_is_E_Field_defined_truth) {
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], _two_D_array[1][two_d_array_pos], _two_D_array[2][two_d_array_pos], threevector(_two_D_array[3][two_d_array_pos], _two_D_array[4][two_d_array_pos], _two_D_array[5][two_d_array_pos], 's'), _is_cartesian, _is_spherical, _is_cylindrical);
                    } else {
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], threevector(_two_D_array[1][two_d_array_pos], _two_D_array[2][two_d_array_pos], _two_D_array[3][two_d_array_pos], 's'));
                    }
                }
                else if (_is_cylindrical) {
		            // Note cylindrical defined in text file as rho, z, phi;
		            // Threevector in order: rho, phi, z
                    if (_is_E_Field_defined_truth) {
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], _two_D_array[1][two_d_array_pos], _two_D_array[2][two_d_array_pos], threevector(_two_D_array[3][two_d_array_pos], _two_D_array[5][two_d_array_pos], _two_D_array[4][two_d_array_pos], 'c'), _is_cartesian, _is_spherical, _is_cylindrical);
                    } else {
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], threevector(_two_D_array[1][two_d_array_pos], _two_D_array[3][two_d_array_pos], _two_D_array[2][two_d_array_pos], 'c'));
                    }
                }
		        k_count+=1;
	        }
	        j_count+=1;
		    k_count=0;
	    }
        i_count+=1;
	    j_count=0;
    }

    _Field_Map_Final = Field_Map_Matrix;

    //Empty because no longer needed
    const std::vector<double> blank = {};
    _two_D_array = {blank};
}

void Field_Map::partition_three_D_grid(){
    const std::vector<std::vector<int>> base_pos_holder= {{0}};
    const std::vector<std::vector<double>> base_value_holder= {{0.0}};
    if (_dimension_val==3){
        find_per_dimension(1);
        find_per_dimension(2);
        find_per_dimension(3);
    } else if (_dimension_val==2){
        find_per_dimension(1);
        find_per_dimension(2);
        _three_d_pos_holder = base_pos_holder;
        _three_d_value_holder = base_value_holder;
    } else if (_dimension_val==1){
        find_per_dimension(1);
        _two_d_pos_holder = base_pos_holder;
        _two_d_value_holder = base_value_holder;
        _three_d_pos_holder = base_pos_holder;
        _three_d_value_holder = base_value_holder;
    } else {
        std::cout<<"WARNING: dimension value unknown."<<std::endl;
    }
}

void Field_Map::find_per_dimension(int dimension){
    int count=1;//initialise with meaningless value since will be replaced
    const std::string warning = "WARNING: invalid dimension provided when seeking to partition grid using find_per_dimension method in Field_Map.cpp";
    if (dimension==1){
        count =_total_one_d_count;
    } else if (dimension==2){
        count = _total_two_d_count;
    } else if (dimension==3){
        count = _total_three_d_count;
    } else{
        std::cout<<warning<<std::endl;
    }
    std::vector<std::vector<double>> value_holder; // holds the value in the field at key positions
    std::vector<std::vector<int>> pos_holder; // holds the index position in the field at key positions
    int num_secs = log2(count); // find the log_2 value
    std::vector<double> line_one_val_holder(3); //line 1
    std::vector<int> line_one_pos_holder = partition(0, count-1); // line 1
    for (int i=0; i<3; i++){
        line_one_val_holder[i] = find_val_from_field_point(line_one_pos_holder[i], dimension);
    }
    pos_holder.push_back(line_one_pos_holder);
    value_holder.push_back(line_one_val_holder);

    // Find the positions in each section from the second line onwards
    for (int i=1; i<num_secs; i++){
	    int num = pos_holder[i-1].size()+ std::pow(2,i);
	    if (num>count/2){
	        // need to go every single 1, this is the last line since a subsequent line would not have unique midpoints.
	        num = count;
	    }
	    std::vector<double> this_iteration_val_holder(num);
	    std::vector<int> this_iteration_pos_holder(num);
	    if (num != count){
	        this_iteration_pos_holder[0] = 0; // first index
	        const int num_to_find= pos_holder[i-1].size()-1;
	        for (int j=0; j<num_to_find; j++){
	            std::vector<int> next_three = partition(pos_holder[i-1][j], pos_holder[i-1][j+1]);
		        this_iteration_pos_holder[2*j+1] = next_three[1];
		        this_iteration_pos_holder[2*j+2] = next_three[2];
	        }
	    } else {
            // this is the last line
	        for (int j=0;j<count;j++){
	            this_iteration_pos_holder[j]=j;
	        }
	    }
        // Fill the values corresponding to the position.
        for (int j=0; j<num; j++){
            this_iteration_val_holder[j] = find_val_from_field_point(this_iteration_pos_holder[j], dimension);
        }
        value_holder.push_back(this_iteration_val_holder);
	    pos_holder.push_back(this_iteration_pos_holder);
    }
    // Update the member variables accordingly
    if(dimension==1){
        _one_d_pos_holder = pos_holder;
	    _one_d_value_holder = value_holder;
    } else if (dimension==2){
        _two_d_pos_holder = pos_holder;
	    _two_d_value_holder = value_holder;
    } else if (dimension==3){
        _three_d_pos_holder = pos_holder;
	    _three_d_value_holder = value_holder;
    } else{
        std::cerr<<warning<<std::endl;
    	throw std::exception();
    }
}

std::vector<int> Field_Map::partition(int start, int end){
    int midway = start+(end-start)/2;
    return {start, midway, end};
}

double Field_Map::find_val_from_field_point(int index, int dimension){
    double value;
    const std::string warning = "WARNING: invalid dimension provided";
    if(dimension==1){
        if (_is_cartesian){
            value = _Field_Map_Final[index][0][0].get_position_x();
        } else if(_is_spherical){
            value = _Field_Map_Final[index][0][0].get_position_r();
        } else if(_is_cylindrical){
            value = _Field_Map_Final[index][0][0].get_position_rho();
        } else {std::cout<<"WARNING: coordinate type not found"<<std::endl;}

    } else if (dimension==2){
        if (_is_cartesian){
            value = _Field_Map_Final[0][index][0].get_position_y();
        } else if(_is_spherical){
            value = _Field_Map_Final[0][index][0].get_position_theta();
        } else if(_is_cylindrical){
            value = _Field_Map_Final[0][index][0].get_position_z();
        } else {std::cout<<"WARNING: coordinate type not found"<<std::endl;}
    } else if (dimension==3){
        if (_is_cartesian){
            value = _Field_Map_Final[0][0][index].get_position_z();
        } else if(_is_spherical){
            value = _Field_Map_Final[0][0][index].get_position_phi();
        } else if(_is_cylindrical){
            value = _Field_Map_Final[0][0][index].get_position_phi();
        } else {std::cout<<"WARNING: coordinate type not found"<<std::endl;}
    } else{
        std::cout<<warning<<std::endl;
    }
    return value;
}

void Field_Map::add_quadratic_fit_details_to_all_points(){
    // The below should already be accounted for
    const std::string error_message = "ERROR: Not enough points provided for dimension to be properly defined. Please provide at least 3 points along each dimension.";
    if (_dimension_val == 1 && _total_one_d_count < 3){std::cout<<error_message;}
    if (_dimension_val == 2  && (_total_one_d_count < 3 || _total_two_d_count<3)){std::cout<<error_message;}
    if (_dimension_val == 3  && (_total_one_d_count < 3 || _total_two_d_count<3 ||_total_three_d_count<3)){std::cout<<error_message;}

    for (int k=0; k<_total_three_d_count; k++){
        for (int j=0; j<_total_two_d_count; j++){
	        for (int i=0; i<_total_one_d_count; i++){
		        if (!_is_E_Field_defined_truth){
		            int pos[3] = {i,j,k};
		            std::vector<std::vector<double>> fits = find_fits(pos);
		            _Field_Map_Final[i][j][k].add_quadratic_fit(fits);
		        }
		        _Field_Map_Final[i][j][k].find_E_field_at_point(_dimension_val, _is_cartesian, _is_spherical, _is_cylindrical, _is_E_Field_defined_truth);
	        }
	    }
    }
}

std::vector<std::vector<double>> Field_Map::find_fits(int position[3]){
    std::vector<std::vector<double>> abc_s;
    // Check dimensionality to ensure correct neighbours are being located
    if (_dimension_val == 1){
        Field_Point points[3];
        int central_pos = position[0];
        // Check if is on the edge and adjust accordingly
        if (central_pos==0){
            central_pos = 1;
        } else if (central_pos==_total_one_d_count-1){
            central_pos = _total_one_d_count-2;
        }
        points[1] =  _Field_Map_Final[central_pos-1][0][0];
        points[0] =  _Field_Map_Final[central_pos][0][0];
        points[2] =  _Field_Map_Final[central_pos+1][0][0];
        double first_D_p0[2];
        double first_D_p1[2];
        double first_D_p2[2];
        first_D_p0[0] = points[0].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical);
    	first_D_p0[1] = points[0].get_value();
        first_D_p1[0] = points[1].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical);
    	first_D_p1[1] = points[1].get_value();
        first_D_p2[0] = points[2].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical);
        first_D_p2[1] = points[2].get_value();
        std::vector<double> first_D_abc = calc_quad_fit(first_D_p0, first_D_p1, first_D_p2);
        abc_s = {first_D_abc, {0.0,0.0,0.0}, {0.0,0.0,0.0}};
    } else if (_dimension_val == 2){
        Field_Point points[6];
        const int one_d_central_pos = position[0];
        const int two_d_central_pos = position[1];
        int shifted_one_d_central_pos = one_d_central_pos;
        int shifted_two_d_central_pos = two_d_central_pos;
        // Check if is on the edge and adjust accordingly
        if (one_d_central_pos==0){
            shifted_one_d_central_pos = 1;
        } else if (one_d_central_pos==_total_one_d_count-1){
            shifted_one_d_central_pos = _total_one_d_count-2;
        }
        if (two_d_central_pos==0){
            shifted_two_d_central_pos = 1;
        } else if (two_d_central_pos==_total_two_d_count-1){
            shifted_two_d_central_pos = _total_two_d_count-2;
        }
        points[0] =  _Field_Map_Final[shifted_one_d_central_pos-1][two_d_central_pos][0];
        points[1] =  _Field_Map_Final[shifted_one_d_central_pos][two_d_central_pos][0];
        points[2] =  _Field_Map_Final[shifted_one_d_central_pos+1][two_d_central_pos][0];
        points[3] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos-1][0];
        points[4] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos][0];
        points[5] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos+1][0];
        // Find position-value pairs (3 pairs in each dimension)
        double first_D_p0[2] = {points[0].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical), points[0].get_value()};
        double first_D_p1[2] = {points[1].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical), points[1].get_value()};
        double first_D_p2[2] = {points[2].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical), points[2].get_value()};
        double second_D_p0[2] = {points[3].get_specified_position(2, _is_cartesian, _is_spherical, _is_cylindrical), points[3].get_value()};
        double second_D_p1[2] = {points[4].get_specified_position(2, _is_cartesian, _is_spherical, _is_cylindrical), points[4].get_value()};
        double second_D_p2[2] = {points[5].get_specified_position(2, _is_cartesian, _is_spherical, _is_cylindrical), points[5].get_value()};
        // Calculate fit through the three position-value pairs in each dimension
        std::vector<double> first_D_abc = calc_quad_fit(first_D_p0, first_D_p1, first_D_p2);
        std::vector<double> second_D_abc = calc_quad_fit(second_D_p0, second_D_p1, second_D_p2);
        abc_s = {first_D_abc, second_D_abc, {0.0,0.0,0.0}};
    } else if (_dimension_val == 3){
        Field_Point points[9];
        const int one_d_central_pos = position[0];
        const int two_d_central_pos = position[1];
        const int three_d_central_pos = position[2];
        int shifted_one_d_central_pos = one_d_central_pos;
        int shifted_two_d_central_pos = two_d_central_pos;
        int shifted_three_d_central_pos = three_d_central_pos;
        // Check if is on the edge and adjust accordingly
        if (one_d_central_pos==0){
            shifted_one_d_central_pos = 1;
        } else if (one_d_central_pos==_total_one_d_count-1){
            shifted_one_d_central_pos = _total_one_d_count-2;
        }
        if (two_d_central_pos==0){
            shifted_two_d_central_pos = 1;
        } else if (two_d_central_pos==_total_two_d_count-1){
            shifted_two_d_central_pos = _total_two_d_count-2;
        }
        if (three_d_central_pos==0){
            shifted_three_d_central_pos = 1;
        } else if (three_d_central_pos==_total_three_d_count-1){
            shifted_three_d_central_pos = _total_three_d_count-2;
        }
        points[0] =  _Field_Map_Final[shifted_one_d_central_pos-1][two_d_central_pos][three_d_central_pos];
        points[1] =  _Field_Map_Final[shifted_one_d_central_pos][two_d_central_pos][three_d_central_pos];
        points[2] =  _Field_Map_Final[shifted_one_d_central_pos+1][two_d_central_pos][three_d_central_pos];
        points[3] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos-1][three_d_central_pos];
        points[4] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos][three_d_central_pos];
        points[5] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos+1][three_d_central_pos];
        points[6] =  _Field_Map_Final[one_d_central_pos][two_d_central_pos][shifted_three_d_central_pos-1];
        points[7] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos][shifted_three_d_central_pos];
        points[8] =  _Field_Map_Final[one_d_central_pos][shifted_two_d_central_pos][shifted_three_d_central_pos+1];
        // Find position-value pairs (3 pairs in each dimension)
        double first_D_p0[2] = {points[0].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical), points[0].get_value()};
        double first_D_p1[2] = {points[1].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical), points[1].get_value()};
        double first_D_p2[2] = {points[2].get_specified_position(1, _is_cartesian, _is_spherical, _is_cylindrical), points[2].get_value()};
        double second_D_p0[2] = {points[3].get_specified_position(2, _is_cartesian, _is_spherical, _is_cylindrical), points[3].get_value()};
        double second_D_p1[2] = {points[4].get_specified_position(2, _is_cartesian, _is_spherical, _is_cylindrical), points[4].get_value()};
        double second_D_p2[2] = {points[5].get_specified_position(2, _is_cartesian, _is_spherical, _is_cylindrical), points[5].get_value()};
        double third_D_p0[2] = {points[6].get_specified_position(3, _is_cartesian, _is_spherical, _is_cylindrical), points[6].get_value()};
        double third_D_p1[2] = {points[7].get_specified_position(3, _is_cartesian, _is_spherical, _is_cylindrical), points[7].get_value()};
        double third_D_p2[2] = {points[8].get_specified_position(3, _is_cartesian, _is_spherical, _is_cylindrical), points[8].get_value()};
        // Calculate fit through the three position-value pairs in each dimension
        std::vector<double> first_D_abc = calc_quad_fit(first_D_p0, first_D_p1, first_D_p2);
        std::vector<double> second_D_abc = calc_quad_fit(second_D_p0, second_D_p1, second_D_p2);
        std::vector<double> third_D_abc = calc_quad_fit(third_D_p0, third_D_p1, third_D_p2);
        abc_s = {first_D_abc, second_D_abc, third_D_abc};
    }
    return abc_s;
}

std::vector<double> Field_Map::calc_quad_fit(double p1[2], double p2[2], double p3[2]){
    // solve in the form B = AC, C is unknown, A is matrix
    // Quadratic of the form y = ax^2+bx +c
    // C = {a,b,c}
    // B = {y1,y2,y3}
    std::vector<double> Y(3);
    // Find differences
    const double allowance = 1e-10;
    const double p1_pos_p2_pos_difference = p1[0]-p2[0];
    const double p2_pos_p3_pos_difference = p2[0]-p3[0];
    const double p1_pos_p3_pos_difference = p1[0]-p3[0];
    const double p1_val_p2_val_difference = p1[1]-p2[1];
    const double p2_val_p3_val_difference = p2[1]-p3[1];
    const double p1_val_p3_val_difference = p1[1]-p3[1];
    const bool p1_pos_equals_p2_pos = std::abs(p1_pos_p2_pos_difference)< allowance;
    const bool p2_pos_equals_p3_pos = std::abs(p2_pos_p3_pos_difference)< allowance;
    const bool p1_pos_equals_p3_pos = std::abs(p1_pos_p3_pos_difference)< allowance;
    const bool p1_val_equals_p2_val = std::abs(p1_val_p2_val_difference)< allowance;
    const bool p2_val_equals_p3_val = std::abs(p2_val_p3_val_difference)< allowance;
    const bool p1_val_equals_p3_val = std::abs(p1_val_p3_val_difference)< allowance;
    if (p1_pos_equals_p2_pos && p2_pos_equals_p3_pos){
	// Solution cannot be found (infinity)
        std::cout<<"ERROR: all three points at same position.\n"
            <<"Fit is incorrectly returned as (0,0,0).\n"
            <<"This occurs at position: "<< p1[0] << std::endl;
        Y = {0.0,0.0,0.0};
    } else if ((p1_pos_equals_p2_pos && !p1_val_equals_p2_val) || (p2_pos_equals_p3_pos && !p2_val_equals_p3_val) || (p1_pos_equals_p3_pos && !p1_val_equals_p3_val) ){
	// Solution cannot be found (infinity)
        std::cout<<"ERROR: two points at same position but not repeated.\n"
            <<"Fit is incorrectly returned as (0,0,0).\n"
            <<"This occurs at position: \n"
            << "\t("<< p1[0]<<", "<< p1[1]<<")\n"
            << "\t("<< p2[0]<<", "<< p2[1]<<")\n"
            << "\t("<< p3[0]<<", "<< p3[1]<<")\n"
            << std::endl;
        Y = {0.0,0.0,0.0};
    } else if ((p1_pos_equals_p2_pos && p1_val_equals_p2_val) || (p1_pos_equals_p3_pos && p1_val_equals_p3_val) ){
	// Repeated point, solution can be found as a straight line
        std::cout<<"Warning: two points are repeated.\n"
            <<"This occurs at position: "<< p1[0] << std::endl;
        const double b = (p2[1]-p3[1])/(p2[0]-p3[0]);
    	const double c = p2[1]-b*p2[0];
    	Y = {0.0,b,c};
    } else if (p2_pos_equals_p3_pos && p2_val_equals_p3_val){
	// Repeated point, solution can be found as a straight line
        std::cout<<"Warning: two points are repeated.\n"
            <<"This occurs at position: "<< p2[0] << std::endl;
        const double b = (p1[1]-p2[1])/(p1[0]-p2[0]);
    	const double c = p1[1]-b*p1[0];
    	Y = {0.0,b,c};
    } else if(p1_val_equals_p2_val && p2_val_equals_p3_val){
	    // This is a horizontal line between all points.
	    Y = {0.0,0.0,p1[1]};
    }
    else if( (p1_val_p2_val_difference)/(p1_pos_p2_pos_difference) == (p2_val_p3_val_difference)/(p2_pos_p3_pos_difference) ){
	// This is a linear fit between all points
	const double b = (p2[1]-p1[1])/(p2[0]-p1[0]);
	const double c = p1[1]-b*p1[0];
	Y = {0.0,b,c};
    }
    else {
    double B[3] = {p1[1], p2[1], p3[1]};
    double A[3][3] = {{p1[0]*p1[0], p1[0], 1}, {p2[0]*p2[0], p2[0], 1}, {p3[0]*p3[0], p3[0], 1}};
    double A_det  = A[0][0]*two_by_two_mat_det(A[1][1],A[1][2],A[2][1],A[2][2])-A[0][1]*two_by_two_mat_det(A[1][0],A[1][2], A[2][0],A[2][2])+A[0][2]*two_by_two_mat_det(A[1][0], A[1][1], A[2][0],A[2][1]);
    // Inverse of matrix A multiplied by A_det
    double top_row[3]={two_by_two_mat_det(A[1][1],A[1][2],A[2][1],A[2][2]),
	               two_by_two_mat_det(A[0][2],A[0][1],A[2][2],A[2][1]),
                       two_by_two_mat_det(A[0][1],A[0][2],A[1][1],A[1][2])};
    double mid_row[3]={two_by_two_mat_det(A[1][2],A[1][0],A[2][2],A[2][0]),
	               two_by_two_mat_det(A[0][0],A[0][2],A[2][0],A[2][2]),
                       two_by_two_mat_det(A[0][2],A[0][0],A[1][2],A[1][0])};
    double bot_row[3]={two_by_two_mat_det(A[1][0],A[1][1],A[2][0],A[2][1]),
	               two_by_two_mat_det(A[0][1],A[0][0],A[2][1],A[2][0]),
                       two_by_two_mat_det(A[0][0],A[0][1],A[1][0],A[1][1])};
    double y1 = (top_row[0]*B[0]+top_row[1]*B[1]+top_row[2]*B[2])/A_det;
    double y2 = (mid_row[0]*B[0]+mid_row[1]*B[1]+mid_row[2]*B[2])/A_det;
    double y3 = (bot_row[0]*B[0]+bot_row[1]*B[1]+bot_row[2]*B[2])/A_det;

    Y = {y1,y2,y3};
    }
    return Y;
}

double Field_Map::two_by_two_mat_det(double upper_left, double upper_right, double bottom_left, double bottom_right){
    return upper_left*bottom_right - upper_right*bottom_left;
}

std::vector<int> Field_Map::find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder, int dimension){
    std::vector<int> closest_positions(2);
    if (value<value_holder[0][0]){
	// out of low bound: state edge position as closest and record out of domain
        closest_positions[0]=0;
	    closest_positions[1]=0;
        _num_times_outside_domain[dimension-1][0]+=1;
    } else if (value>value_holder[0][2]){
	//out of high bound: state edge position as closest and record out of domain
        closest_positions[0]=pos_holder[0][2];
	    closest_positions[1]=pos_holder[0][2];
	    _num_times_outside_domain[dimension-1][1]+=1;
    } else{
        int test_point = 1;
        int num_secs = pos_holder.size();
        for (int i=0; i<num_secs-2; i++){
	        if (value > value_holder[i][test_point]){
                test_point=2*test_point+1; // for the next line, the test point is shifted
	        } else {
	            test_point=2*test_point-1;
	        }
        }
        // Penultimate line:
        //     check partition above or below depending on
        if (value < value_holder[num_secs-2][test_point]){
            test_point -= 1;
        }
        // check diff to this value, diff to next value etc. until buffer
        int num_to_check = 1+pos_holder[num_secs-2][test_point+1] - pos_holder[num_secs-2][test_point];
	    int index_of_test_point = pos_holder[num_secs-2][test_point];
        int closest_pos = index_of_test_point;
	    double closest_diff = value-value_holder[num_secs-2][test_point];
	    int second_closest_pos = index_of_test_point+1;
	    for (int i=1; i<num_to_check;i++){
	        double this_diff = value-value_holder[num_secs-1][index_of_test_point+i];
	        if (std::abs(this_diff)<std::abs(closest_diff)){
	            if(this_diff>0){
		            second_closest_pos = closest_pos+2;
		        } else {
		            second_closest_pos = closest_pos;
		        }
		        closest_pos = index_of_test_point + i;
		        closest_diff = this_diff;
	        }
	    }
        closest_positions[0] = closest_pos;
	    closest_positions[1] = second_closest_pos;
    }
    return closest_positions;
}

threevector Field_Map::find_approx_value(threevector point){
    // Begin by finding the close points in each dimension
    // Store the position indexes of the near points
    std::vector<int> one_d_close_pos(2);
    std::vector<int> two_d_close_pos(2);
    std::vector<int> three_d_close_pos(2);
    // number of points for 1 dimension = 2; 2d=4;3d=8
    const int num_near_points= std::pow(2, _dimension_val);
    // Store the actual field points of the near points
    Field_Point near_points[num_near_points];
    if(_dimension_val==1){
         if (_is_cartesian){
             one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder, 1);
         } else if (_is_spherical){
	         one_d_close_pos = find_closest(point.mag3(), _one_d_pos_holder, _one_d_value_holder, 1);
         } else if (_is_cylindrical){
             one_d_close_pos = find_closest(point.getrho(), _one_d_pos_holder, _one_d_value_holder, 1);
         }
         near_points[0] = _Field_Map_Final[one_d_close_pos[0]][0][0];
         near_points[1] = _Field_Map_Final[one_d_close_pos[1]][0][0];
    } else if (_dimension_val==2){
        if (_is_cartesian){
            one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder, 1);
            two_d_close_pos = find_closest(point.gety(), _two_d_pos_holder, _two_d_value_holder, 2);
        } else if (_is_spherical){
            one_d_close_pos = find_closest(point.mag3(), _one_d_pos_holder, _one_d_value_holder, 1);
            two_d_close_pos = find_closest(point.gettheta(), _two_d_pos_holder, _two_d_value_holder, 2);
        } else if (_is_cylindrical){
            one_d_close_pos = find_closest(point.getrho(), _one_d_pos_holder, _one_d_value_holder, 1);
            two_d_close_pos = find_closest(point.getz(), _two_d_pos_holder, _two_d_value_holder, 2);
        }
        near_points[0] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][0];
        near_points[1] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[0]][0];
        near_points[2] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[1]][0];
        near_points[3] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[1]][0];
    } else if (_dimension_val==3){
        if (_is_cartesian){
            one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder, 1);
            two_d_close_pos = find_closest(point.gety(), _two_d_pos_holder, _two_d_value_holder, 2);
            three_d_close_pos = find_closest(point.getz(), _three_d_pos_holder, _three_d_value_holder, 3);
        } else if (_is_spherical){
            one_d_close_pos = find_closest(point.mag3(), _one_d_pos_holder, _one_d_value_holder, 1);
            two_d_close_pos = find_closest(point.gettheta(), _two_d_pos_holder, _two_d_value_holder, 2);
            three_d_close_pos = find_closest(point.getphi(), _three_d_pos_holder, _three_d_value_holder, 3);
        } else if (_is_cylindrical){
            one_d_close_pos = find_closest(point.getrho(), _one_d_pos_holder, _one_d_value_holder, 1);
            two_d_close_pos = find_closest(point.getz(), _two_d_pos_holder, _two_d_value_holder, 2);
            three_d_close_pos = find_closest(point.getphi(), _three_d_pos_holder, _three_d_value_holder, 3);
        }
        near_points[0] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][three_d_close_pos[0]];
        near_points[1] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[0]][three_d_close_pos[0]];
        near_points[2] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[1]][three_d_close_pos[0]];
        near_points[3] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[1]][three_d_close_pos[0]];
        near_points[4] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][three_d_close_pos[1]];
        near_points[5] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[0]][three_d_close_pos[1]];
        near_points[6] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[1]][three_d_close_pos[1]];
        near_points[7] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[1]][three_d_close_pos[1]];
    } else {std::cout<<"WARNING: unknown dimension number"<<std::endl;}

    // Now the closest field points are found, find and average the electric field at each point
    double summed_weight = 0;
    threevector E_field_values[num_near_points];
    double weights[num_near_points];
    threevector summed_weighted_values = threevector(0,0,0);
    bool is_exact_pos_found = false;
    threevector exact_E_field_val;
    for (int i=0;i<num_near_points;i++){
        if (_is_E_Field_defined_truth || !_CALCULATE_NEARBY_E_FIELD_ON_CURVES) {
            E_field_values[i] = near_points[i].get_electric_field_at_point();
        } else {
            E_field_values[i] = near_points[i].get_nearby_E_field(point, _dimension_val, _is_cartesian, _is_spherical, _is_cylindrical, _is_E_Field_defined_truth);
        }
        // For averaging, see how far the nearby field point is from the desired
        //     true position
        double distance_squared;
        if (_dimension_val==1){
            if (_is_cartesian){
                distance_squared = std::pow(near_points[i].get_position_x()-point.getx(),2);
            } else if (_is_spherical) {
                distance_squared = std::pow(near_points[i].get_position_r()-point.mag3(),2);
            } else if (_is_cylindrical){
                distance_squared = std::pow(near_points[i].get_position_rho()-point.getrho(),2);
            }
        } else if (_dimension_val==2){
            if (_is_cartesian){
                distance_squared = std::pow(near_points[i].get_position_x()-point.getx(),2)+std::pow(near_points[i].get_position_y()-point.gety(),2);
            } else if (_is_spherical) {
                distance_squared = std::pow(near_points[i].get_position_r()-point.mag3(),2) + std::pow(near_points[i].get_position_theta()-point.gettheta(),2);
            } else if (_is_cylindrical){
                distance_squared = std::pow(near_points[i].get_position_rho()-point.getrho(),2) + std::pow(near_points[i].get_position_z()-point.getz(),2);
            }
        } else if (_dimension_val==3){
            if (_is_cartesian){
                distance_squared = std::pow(near_points[i].get_position_x()-point.getx(),2)+std::pow(near_points[i].get_position_y()-point.gety(),2)+std::pow(near_points[i].get_position_z()-point.getz(),2);
            } else if (_is_spherical) {
                distance_squared = std::pow(near_points[i].get_position_r()-point.mag3(),2) + std::pow(near_points[i].get_position_theta()-point.gettheta(),2) + std::pow(near_points[i].get_position_phi()-point.getphi(),2);
            } else if (_is_cylindrical){
                distance_squared = std::pow(near_points[i].get_position_rho()-point.getrho(),2) + std::pow(near_points[i].get_position_z()-point.getz(),2) + std::pow(near_points[i].get_position_phi()-point.getphi(),2);
            }
        }

	    if (distance_squared == 0){
	        is_exact_pos_found = true;
	        exact_E_field_val = E_field_values[i];
	    }
	    else {
	        weights[i] = 1/distance_squared;
	        summed_weight += weights[i];
	        summed_weighted_values += E_field_values[i]*weights[i];
	    }
    }
    threevector approximate_E_Field_val;
    if (is_exact_pos_found){
        approximate_E_Field_val = exact_E_field_val;
    }
    else {
        approximate_E_Field_val = threevector(summed_weighted_values.getx()/summed_weight, summed_weighted_values.gety()/summed_weight, summed_weighted_values.getz()/summed_weight);
    }
    return approximate_E_Field_val;
}


void Field_Map::check_coordinate_domain_symmetry(int dimension_to_check){
    double min;
    double max;
    const double delta = 1e-4; // gives user input a little margin for error
    const std::string coordinate_name = get_coordinate_name(dimension_to_check);
    const std::string symmetric_warning_message = "Warning! The user-defined" + coordinate_name + "range is not symmetric. Proceeding, but take caution.";
    if (dimension_to_check == 1){
         min = _one_d_value_holder[0][0];
         max = _one_d_value_holder[0][2];
    } else if (dimension_to_check == 2){
         min = _two_d_value_holder[0][0];
         max = _two_d_value_holder[0][2];
    } else if (dimension_to_check == 3){
         min = _three_d_value_holder[0][0];
         max = _three_d_value_holder[0][2];
    }
    if (std::abs(min+max)>delta){
	    std::cout<<symmetric_warning_message<<std::endl;
	}
}

void Field_Map::check_coordinate_domain_completeness(int dimension_to_check){
    double min;
    double max;
    const double delta = 1e-4; // gives user input a little margin for error
    const std::string coordinate_name = get_coordinate_name(dimension_to_check);
    const std::string symmetric_warning_message = "\nWarning! The user-defined" + coordinate_name + "range is not complete. Proceeding using extrapolation, but take caution.";
    if (dimension_to_check == 2){
        min = _two_d_value_holder[0][0];
        max = _two_d_value_holder[0][2];
        if (min > delta){
            std::cout<<symmetric_warning_message<<"\n\tDefined lower limit: "<<min<<std::endl;
	    }
	    if (max < PI-delta){
            std::cout<<symmetric_warning_message<<"\n\tDefined upper limit: "<<max<<std::endl;
	    }
    } else if (dimension_to_check == 3){
        min = _three_d_value_holder[0][0];
        max = _three_d_value_holder[0][2];
        if (min > delta){
            std::cout<<symmetric_warning_message<<"\n\tDefined lower limit: "<<min<<std::endl;
	    }
	if (max < 2*PI-delta){
            std::cout<<symmetric_warning_message<<"\n\tDefined upper limit: "<<max<<std::endl;
	    }
    }
}

void Field_Map::find_domain_limits(){
    const double MAX = 1e6; // Large value, used for situation where dimension is not defined by custom map.
    if (_is_cartesian){
	    const double limiting_x = std::min(std::abs(_one_d_value_holder[0][0]), _one_d_value_holder[0][2]);
	    if (_dimension_val==1){
	        _map_z_min = -MAX;
            _map_z_max = MAX;
	        _map_rho_max = limiting_x;
	    } else {
	        const double limiting_y = std::min(std::abs(_two_d_value_holder[0][0]), _two_d_value_holder[0][2]);
	        _map_rho_max = std::sqrt(limiting_x*limiting_x+limiting_y*limiting_y);
	        if (_dimension_val==2){
	            _map_z_min = -MAX;
                _map_z_max = MAX;
	        } else if (_dimension_val==3){
		        _map_z_min = _three_d_value_holder[0][0];
                _map_z_max = _three_d_value_holder[0][2];
	        }
	    }
    } else if (_is_spherical) {
        const double map_r = _one_d_value_holder[0][2];
        if (_is_spherical_injection){
            _map_z_max = map_r;
            _map_z_min = -map_r;
            _map_rho_max = map_r;
        }
        if (_is_cylindrical_injection){
	       // Find the largest cylinder that will fit within sphere
	        _map_z_max = map_r/std::sqrt(3);
	        _map_z_min = -_map_z_max;
            _map_rho_max = _map_z_max*std::sqrt(2);
        }
    } else if (_is_cylindrical){
        _map_rho_max = _one_d_value_holder[0][2];
	    if (_dimension_val==1){
	        _map_z_min = -MAX;
            _map_z_max = MAX;
	    } else {
	        _map_z_min = _two_d_value_holder[0][0];
            _map_z_max = _two_d_value_holder[0][2];
	    }
    }
}


void Field_Map::set_new_map_limits(double new_rho_limit){
    if (_is_spherical && _is_cylindrical_injection){
	    // Find the largest cylinder that will fit within sphere
        const double buffer = 1.1;
	    const double map_r = _one_d_value_holder[0][2];
	    _map_z_max = std::sqrt(map_r*map_r - new_rho_limit*new_rho_limit)/buffer;
	    _map_z_min = -_map_z_max;
	    _map_rho_max = new_rho_limit;
    }
}

void Field_Map::check_dust_grain_map_coverage(){
    // Just check spherical and cylindrical since Cartesian should be fine, so long as it is symmetric and covers dimpl domain which is checked elsewhere.
    if (_is_spherical || _is_cylindrical) {
        const double limiting_a = std::min(std::min(_a1, _a2), _a3);
        const double map_r_min = _one_d_value_holder[0][0];
        if (limiting_a < map_r_min){
            const std::string limiting_a_string = std::to_string(limiting_a);
            const std::string map_r_min_string = std::to_string(map_r_min);
            const std::string r_limit_warning_message = "Warning! The user-defined custom potential map does not reach the dust grain surface. Proceeding using extrapolation, but please take caution.";

            std::cout<<r_limit_warning_message<<"\n\t Dust grain radius limiting value: "<<limiting_a_string<<"\n\t Lowest radial potential value: "<<map_r_min_string <<std::endl;
        }
    }
}

std::string Field_Map::get_coordinate_name(int dimension){
    std::string coordinate_name;
    if (_is_cartesian){
        if (dimension==1){coordinate_name = "x";}
        else if (dimension==2){coordinate_name = "y";}
        else if (dimension==3){coordinate_name = "z";}
    } else if (_is_spherical){
        if (dimension==1){coordinate_name = "r";}
        else if (dimension==2){coordinate_name = "theta";}
        else if (dimension==3){coordinate_name = "phi";}
    } else if (_is_cylindrical){
        if (dimension==1){coordinate_name = "rho";}
        else if (dimension==2){coordinate_name = "z";}
        else if (dimension==3){coordinate_name = "phi";}
    }
    return coordinate_name;
}

void Field_Map::check_domain(){
   // Check for symmetry and completeness.
   // Check for logical choice of coordinate system.
   // Check the z and radial injection limits set by dimpl are covered by the custom map
   // Check for dust limits
   const std::string conflicting_map_injection_warning_message = "Warning! The injection method coordinate system chosen is different to the custom potential map coordinate system.\n\tProceeding, but please be mindful.";

   if (_is_cartesian){
       if (_is_spherical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tSpherical injection, cartesian potential map."<<std::endl;
       }
       if (_is_cylindrical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tCylindrical injection, cartesian potential map."<<std::endl;
       }
       check_coordinate_domain_symmetry(1);
       if (_dimension_val==2||_dimension_val==3){
           check_coordinate_domain_symmetry(2);
           if (_dimension_val==3){
               check_coordinate_domain_symmetry(3);
           }
       }
   } else if (_is_spherical){
       if (_is_cylindrical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tCylindrical injection, spherical potential map."<<std::endl;
       }
       if (_dimension_val==2||_dimension_val==3){
           check_coordinate_domain_completeness(2);
           if (_dimension_val==3){
               check_coordinate_domain_completeness(3);
           }
       }
   } else if (_is_cylindrical){
       if (_is_spherical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tSpherical injection, cylindrical potential map."<<std::endl;
       }
       if (_dimension_val==2||_dimension_val==3){
           check_coordinate_domain_symmetry(2);
           if (_dimension_val==3){
               check_coordinate_domain_completeness(3);
           }
       }
   }
   check_dust_grain_map_coverage();
}


int Field_Map::get_shrunk_point(double new_val, double old_val, int dimension, bool are_values_positive){
    int new_pos;
    // Find the appropriate position and value holders for this dimension
    std::vector<std::vector<int>> pos_holder;
    std::vector<std::vector<double>> value_holder;
    if (dimension == 1){
        pos_holder = _one_d_pos_holder;
	    value_holder = _one_d_value_holder;
    } else if (dimension == 2){
        pos_holder = _two_d_pos_holder;
	    value_holder = _two_d_value_holder;
    } else if (dimension == 3){
        pos_holder = _three_d_pos_holder;
	    value_holder = _three_d_value_holder;
    }
    if (are_values_positive){
        if (new_val < old_val){
	        std::vector<int> close_pos= find_closest(new_val, pos_holder, value_holder, dimension);
	        new_pos = std::max(close_pos[0], close_pos[1]);
	    } else {
	        new_pos = pos_holder[0][2];
	    }
    } else {
        if (new_val > old_val){
	        std::vector<int> close_pos= find_closest(new_val, pos_holder, value_holder, dimension);
	        new_pos = std::min(close_pos[0], close_pos[1]);
	    } else {
	        new_pos = pos_holder[0][0];
	    }
    }
    return new_pos;
}

bool Field_Map::check_if_vectors_differ(std::vector<int> one, std::vector<int> two){
    bool is_different = false;
    while (!is_different){
        int vec_one_size = one.size();
	for (int i=0; i<vec_one_size; i++){
	        if (one[i]!=two[i]){is_different = true;}
	    }
    }
    return is_different;
}

void Field_Map::shrink_map_to_fit(double z_min, double z_max, double rho_max){
    const int buffer = 2; // add a buffer to ensure that the new domain fully fits within the shrunk map
    // assign with old
    int new_min_dim1_pos = _one_d_pos_holder[0][0];
    int new_max_dim1_pos = _one_d_pos_holder[0][2];
    int new_min_dim2_pos;
    int new_max_dim2_pos;
    int new_min_dim3_pos;
    int new_max_dim3_pos;
    if (_dimension_val==1){
        new_min_dim2_pos = 0;
        new_max_dim2_pos = 0;
        new_min_dim3_pos = 0;
        new_max_dim3_pos = 0;
    } else {
        new_min_dim2_pos = _two_d_pos_holder[0][0];
        new_max_dim2_pos = _two_d_pos_holder[0][2];
	    if (_dimension_val==2){
	        new_min_dim3_pos = 0;
            new_max_dim3_pos = 0;
	    } else {
	        new_min_dim3_pos = _three_d_pos_holder[0][0];
            new_max_dim3_pos = _three_d_pos_holder[0][2];
	    }
    }
    bool is_shrinking_needed;
    // Find new points at which the limits would be encapsulated
    if(_is_cartesian){
	    new_min_dim1_pos = get_shrunk_point(-rho_max, _one_d_value_holder[0][0], 1, false);
	    new_max_dim1_pos = get_shrunk_point(rho_max, _one_d_value_holder[0][2], 1, true);
	    if (_dimension_val==2 || _dimension_val==3){
	        new_min_dim2_pos = get_shrunk_point(-rho_max, _two_d_value_holder[0][0], 2, false);
	        new_max_dim2_pos = get_shrunk_point(rho_max, _two_d_value_holder[0][2], 2, true);
	        if (_dimension_val==3){
		        new_min_dim3_pos = get_shrunk_point(z_min, _map_z_min, 3, false);
	            new_max_dim3_pos = get_shrunk_point(z_max, _map_z_max, 3, true);
	        }
	    }
    } else if(_is_spherical){
	    double limiting_new_val;
	    const double limiting_z = std::max(-z_min, z_max);
	    if (_is_spherical_injection){
	        limiting_new_val = std::max(limiting_z, rho_max);
	    } else if (_is_cylindrical_injection){
	        limiting_new_val = std::sqrt(limiting_z*limiting_z + rho_max*rho_max);
	    }
	    new_max_dim1_pos = get_shrunk_point(limiting_new_val, _one_d_value_holder[0][2], 1, true);
    } else if(_is_cylindrical){
	    new_max_dim1_pos = get_shrunk_point(rho_max, _one_d_value_holder[0][2], 1, true);
	    if (_dimension_val==2 || _dimension_val==3){
	        new_min_dim2_pos = get_shrunk_point(z_min, _two_d_value_holder[0][0], 2, false);
	        new_max_dim2_pos = get_shrunk_point(z_max, _two_d_value_holder[0][2], 2, true);
	    }
    }
    // Check if the new positions differ to the old positions and state shrinking required if they are
    if (_dimension_val==1){
	    std::vector<int> old_positions = {_one_d_pos_holder[0][0], _one_d_pos_holder[0][2]};
	    std::vector<int> new_positions = {new_min_dim1_pos, new_max_dim1_pos};
	    is_shrinking_needed = check_if_vectors_differ(old_positions,new_positions);
    } else if (_dimension_val==2){
	    std::vector<int> old_positions = {_one_d_pos_holder[0][0], _one_d_pos_holder[0][2], _two_d_pos_holder[0][0], _two_d_pos_holder[0][2]};
	    std::vector<int> new_positions = {new_min_dim1_pos, new_max_dim1_pos, new_min_dim2_pos, new_max_dim2_pos};
	    is_shrinking_needed = check_if_vectors_differ(old_positions,new_positions);
    } else if (_dimension_val==3){
	    std::vector<int> old_positions = {_one_d_pos_holder[0][0], _one_d_pos_holder[0][2], _two_d_pos_holder[0][0], _two_d_pos_holder[0][2], _three_d_pos_holder[0][0], _three_d_pos_holder[0][2]};
	    std::vector<int> new_positions = {new_min_dim1_pos, new_max_dim1_pos, new_min_dim2_pos, new_max_dim2_pos, new_min_dim3_pos, new_max_dim3_pos};
	    is_shrinking_needed = check_if_vectors_differ(old_positions,new_positions);
    }
    if (is_shrinking_needed){
        // Find the new counts
	    int one_d_count = std::min(buffer + new_max_dim1_pos - new_min_dim1_pos,_total_one_d_count);
	    int two_d_count = std::min(buffer + new_max_dim2_pos - new_min_dim2_pos,_total_two_d_count);
	    int three_d_count = std::min(buffer + new_max_dim3_pos - new_min_dim3_pos,_total_three_d_count);
        // if the buffer takes it past the count, then no need to shrink
	// NOTE: below is not working
	    if(!(one_d_count==_total_one_d_count && two_d_count==_total_two_d_count && three_d_count==_total_three_d_count)){
	        std::vector<std::vector<std::vector<Field_Point>>> New_Field_Map_Matrix(one_d_count, std::vector<std::vector<Field_Point>>(two_d_count, std::vector<Field_Point>(three_d_count)));
	        int i_count = 0;
	        int j_count = 0;
	        int k_count = 0;
	        for (int i=new_min_dim1_pos; i<new_max_dim1_pos+1; i++){
	            for (int j=new_min_dim2_pos; j<new_max_dim2_pos+1; j++){
		            for (int k=new_min_dim3_pos; k<new_max_dim3_pos+1; k++){
		                New_Field_Map_Matrix[i_count][j_count][k_count] = _Field_Map_Final[i][j][k];
		                k_count += 1;
		            }
		            j_count += 1;
		            k_count = 0;
	            }
	            i_count += 1;
	            j_count = 0;
	            k_count = 0;
	        }
	        // Update accordingly
            _total_one_d_count = one_d_count;
            _total_two_d_count = two_d_count;
            _total_three_d_count = three_d_count;
	        _Field_Map_Final = New_Field_Map_Matrix;
	        partition_three_D_grid();
            check_domain();
            find_domain_limits();
        }
    }
}

void Field_Map::check_num_times_outside_domain(){
    const std::string domain_breech_warning = "Warning! Instance of outside of custom-defined domain found.\n\t";
    if (_num_times_outside_domain[0][0]!=0){
	std::cout<<domain_breech_warning<<_num_times_outside_domain[0][0]
	<<" instances found below the First Dimension range."
	<<std::endl;
    } if (_num_times_outside_domain[0][1]!=0){
	std::cout<<domain_breech_warning<<_num_times_outside_domain[0][0]
	<<" instances found above the First Dimension range."
	<<std::endl;
    } if (_num_times_outside_domain[1][0]!=0){
	std::cout<<domain_breech_warning<<_num_times_outside_domain[0][0]
	<<" instances found below the Second Dimension range."
	<<std::endl;
    } if (_num_times_outside_domain[1][1]!=0){
	std::cout<<domain_breech_warning<<_num_times_outside_domain[0][0]
	<<" instances found above the Second Dimension range."
	<<std::endl;
    } if (_num_times_outside_domain[2][0]!=0){
	std::cout<<domain_breech_warning<<_num_times_outside_domain[0][0]
	<<" instances found below the Third Dimension range."
	<<std::endl;
    } if (_num_times_outside_domain[2][1]!=0){
	std::cout<<domain_breech_warning<<_num_times_outside_domain[0][0]
	<<" instances found above the Third Dimension range."
	<<std::endl;
    }
}
