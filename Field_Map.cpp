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
Field_Map::Field_Map(std::string field_file_name)
{
    std::clock_t start;
    double duration;
    start = std::clock();
    std::cout<<"Opened constructor"<<std::endl;
    text_to_string_vector(field_file_name);
    duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
    std::cout<<"text to string vector complete, time: "<< duration <<std::endl;
    find_header_details();
    duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
    std::cout<<"Header details found, time: "<<duration<<std::endl;
    convert_data_line_to_multi_D_array();
    duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
    std::cout<<"Convereted data to multi d array, time: "<<duration<<std::endl;
    partition_three_D_grid();
    duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
    std::cout<<"Partitioned three D grid, time: "<< duration<<std::endl;
    if (!_is_E_Field_defined_truth){
        add_quadratic_fit_details_to_all_points();
        duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
        std::cout<<"Fitted all points with a quadratic fit, time: "<<duration<<std::endl;
    }
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
    _ordered_truth = string_to_bool(split_after_delim(_ordered_delim, _line_vector[3])); // and convert string to boolean
    _is_E_Field_defined_truth = string_to_bool(split_after_delim(_is_E_Field_defined_delim, _line_vector[4])); // and convert string to boolean
    // Check the overriding coordinate system.
    const std::string dimension_mismatch_warning = "WARNING: coordinate system dimensions do not match those specified under `Number of Dimensions'. Proceeding but be mindful of potential errors.";
    if (_coord_type == _CART_1 || _coord_type == _CART_2 || _coord_type == _CART_3){
        _is_cartesian = true;
        if (_coord_type == _CART_1 && _dimension_val!=1){
            std::cout<<dimension_mismatch_warning<<std::endl;
        } else if (_coord_type == _CART_2 && _dimension_val!=2){
            std::cout<<dimension_mismatch_warning<<std::endl;
        } else if (_coord_type == _CART_3 && _dimension_val!=3){
            std::cout<<dimension_mismatch_warning<<std::endl;
        }
    } else if (_coord_type == _SPHER_1 || _coord_type == _SPHER_2 || _coord_type == _SPHER_3){
        _is_spherical = true;
        if (_coord_type == _SPHER_1 && _dimension_val!=1){
            std::cout<<dimension_mismatch_warning<<std::endl;
        } else if (_coord_type == _SPHER_2 && _dimension_val!=2){
            std::cout<<dimension_mismatch_warning<<std::endl;
        } else if (_coord_type == _SPHER_3 && _dimension_val!=3){
            std::cout<<dimension_mismatch_warning<<std::endl;
        }
    } else if (_coord_type == _CYL_1 || _coord_type == _CYL_2 || _coord_type == _CYL_3){
        _is_cylindrical = true;
        if (_coord_type == _CYL_1 && _dimension_val!=1){
            std::cout<<dimension_mismatch_warning<<std::endl;
        } else if (_coord_type == _CYL_2 && _dimension_val!=2){
            std::cout<<dimension_mismatch_warning<<std::endl;
        } else if (_coord_type == _CYL_3 && _dimension_val!=3){
            std::cout<<dimension_mismatch_warning<<std::endl;
        }
    } else {
        std::cout<<"WARNING: coordinate system is not recognised, assumed Cartesian."<<std::endl;
        _is_cartesian = true;
    }
}

std::string Field_Map::split_after_delim(std::string delim, std::string str_to_split)
{
    const int end_delim_pos = str_to_split.find(_end_delim);
    const int start_index = str_to_split.find(delim) + delim.length();
    return str_to_split.substr(str_to_split.find(delim) + delim.length(), end_delim_pos-start_index);
}

bool Field_Map::string_to_bool(std::string string){
    const std::string unknown_boolean_warning = "WARNING: ordered truth of custom field file is unknown, assuming ordered.";
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
	        std::cout<< unknown_boolean_warning <<std::endl;
	        string_as_boolean = true;
	    }
    } else {
        // the text boolean case
        std::transform(string.begin(), string.end(), string.begin(), [](unsigned char c) {return std::tolower(c); }); // convert string to lower case
	    if (string=="true"){string_as_boolean = true;}
	    else if (string=="false"){string_as_boolean = false;}
        else {
	        std::cout << unknown_boolean_warning << std::endl;
		std::cout<<"\t String stated to be: "<<string<<std::endl;
	        string_as_boolean = true;
	    }
    }
    return string_as_boolean;
}


void Field_Map::convert_data_line_to_multi_D_array(){
    // Begin by finding the size
    double two_d_array_width;
    if (_is_E_Field_defined_truth) {_num_data_points_per_line=2*_dimension_val; two_d_array_width=6;}
    else {_num_data_points_per_line=1+_dimension_val;two_d_array_width=4;}

    _num_data_lines = _line_vector.size() - _num_lines_in_header;

    _two_D_array = std::vector<std::vector<double>>(two_d_array_width, std::vector<double>(_num_data_lines));
    for (int i=0; i<_num_data_lines; i++){
	    std::string line = _line_vector[i+_num_lines_in_header];
	    std::string line_data[_num_data_points_per_line];
	    double line_values[_num_data_points_per_line];
	    for (int j=0; j<_num_data_points_per_line; j++){
	        int pos = line.find(_data_delim);
	        line_data[j] = line.substr(0, pos);
	        line.erase(0, pos + _data_delim.length());
		line_values[j] = 0.0 + atof(line_data[j].c_str());
	    }
        if (_dimension_val==1){
            if (_is_E_Field_defined_truth){
		_two_D_array[0][i] = line_values[0];
                _two_D_array[1][i] = 0.0;
                _two_D_array[2][i] = 0.0;
		_two_D_array[3][i] = line_values[1];
                _two_D_array[4][i] = 0.0;
                _two_D_array[5][i] = 0.0;
            }
            else {
		_two_D_array[0][i] = line_values[0];
		_two_D_array[1][i] = line_values[1];
                _two_D_array[2][i] = 0.0;
                _two_D_array[3][i] = 0.0;
            }
        }
        if (_dimension_val==2){
            if (_is_E_Field_defined_truth){
		_two_D_array[0][i] = line_values[0];
		_two_D_array[1][i] = line_values[1];
                _two_D_array[2][i] = 0.0;
		_two_D_array[3][i] = line_values[2];
		_two_D_array[4][i] = line_values[3];
                _two_D_array[5][i] = 0.0;
            }
            else {
		_two_D_array[0][i] = line_values[0];
		_two_D_array[1][i] = line_values[1];
		_two_D_array[2][i] = line_values[2];
                _two_D_array[3][i] = 0.0;
            }
        }
        if (_dimension_val==3){
            if (_is_E_Field_defined_truth){
		_two_D_array[0][i] = line_values[0];
		_two_D_array[1][i] = line_values[1];
                _two_D_array[2][i] = line_values[2];
		_two_D_array[3][i] = line_values[3];
		_two_D_array[4][i] = line_values[4];
                _two_D_array[5][i] = line_values[5];
            }
            else {
		_two_D_array[0][i] = line_values[0];
		_two_D_array[1][i] = line_values[1];
		_two_D_array[2][i] = line_values[2];
                _two_D_array[3][i] = line_values[3];
            }
        }
    }
    // Convert grid representation to multi dimensional array of Field Points
    // Note to self: only working for three dimensional
    //
    // Begin by finding number of points in each dimension
    int one_d_count = 1;
    int two_d_count = 1;
    int three_d_count = 1;
    int first_dimension_pos = 1;
    if (_is_E_Field_defined_truth){first_dimension_pos = 3;}
    double memory_holder = _two_D_array[first_dimension_pos][0]; // start off with the first value in the first dimension
    if (!_ordered_truth){std::cout<<"WARNING: method only currently works for ordered files. Proceeding, but be mindful of potential errors."<<std::endl;}
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
                if (third_d_current_val<memory_holder+1.0e-300){// note occasionally a rounding error hence this tiny shift
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
		if (two_d_pos < 0.0){
		    two_d_low_error_count+=1;
		} else if (two_d_pos >PI){
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


    // Populate multi dimensional array
    std::vector<std::vector<std::vector<Field_Point>>> Field_Map_Matrix(one_d_count, std::vector<std::vector<Field_Point>>(two_d_count, std::vector<Field_Point>(three_d_count)));
    int i_count = 0;
    int j_count = 0;
    int k_count = 0;
    const double delta = 1e-15;
    bool is_theta_zero_found = false;
    //threevector E_Field;
    //threevector pos;
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
			// Check if theta is 0 (in this case phi would be badly defined so shift theta by delta)
			double E_field_theta = _two_D_array[1][two_d_array_pos];
			double pos_theta = _two_D_array[4][two_d_array_pos];
		        if ( E_field_theta == 0.0) {is_theta_zero_found = true; E_field_theta+=delta;}
		        if ( pos_theta == 0.0) {is_theta_zero_found = true; pos_theta+=delta;}

                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], E_field_theta, _two_D_array[2][two_d_array_pos], threevector(_two_D_array[3][two_d_array_pos], pos_theta, _two_D_array[5][two_d_array_pos], 's'), _is_cartesian, _is_spherical, _is_cylindrical);
                    } else {
			double pos_theta = _two_D_array[2][two_d_array_pos];
		        if ( pos_theta == 0.0) {is_theta_zero_found = true; pos_theta+=delta;}
                        Field_Map_Matrix[i_count][j_count][k_count] = Field_Point(_two_D_array[0][two_d_array_pos], threevector(_two_D_array[1][two_d_array_pos], pos_theta, _two_D_array[3][two_d_array_pos], 's'));
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
    if (is_theta_zero_found){std::cout<<"WARNING: A provided theta value was 0 which would cause a badly defined phi.\n\t Theta has been shifted by "<<delta<<" to accomodate this, but caution should be taken."<<std::endl;}

    _Field_Map_Final = Field_Map_Matrix;

    //Empty because no longer needed
    const std::vector<double> blank = {};
    _two_D_array = {blank};
}

void Field_Map::partition_three_D_grid(){
    // Note to self: check how many dimensions are required
    // vector of vectors, in one direction the partition values, in the other direction additional granularity
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
    int count=1;//initialise with meaningless value
    const std::string warning = "WARNING: invalid dimension provided";
    if(dimension==1){
        count =_total_one_d_count;
    }else if (dimension==2){
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
        std::cout<<warning<<std::endl;
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
    // Check properly defined
    const std::string error_message = "ERROR: Not enough points provided for dimension to be properly defined. Please provide at least 3 points along each dimension.";
    if (_dimension_val == 1 && _total_one_d_count < 3){std::cout<<error_message;}
    if (_dimension_val == 2  && (_total_one_d_count < 3 || _total_two_d_count<3)){std::cout<<error_message;}
    if (_dimension_val == 3  && (_total_one_d_count < 3 || _total_two_d_count<3 ||_total_three_d_count<3)){std::cout<<error_message;}
    for (int k=0; k<_total_three_d_count; k++){
        for (int j=0; j<_total_two_d_count; j++){
	        for (int i=0; i<_total_one_d_count; i++){
		        int pos[3] = {i,j,k};
		        std::vector<std::vector<double>> fits = find_fits(pos);
		        _Field_Map_Final[i][j][k].add_quadratic_fit(fits);
		        // Note to self: below needed?
			//_Field_Map_Final[i][j][k].find_E_field_at_point(_dimension_val, _is_cartesian, _is_spherical, _is_cylindrical);
	        }
	    }
    }
}

std::vector<std::vector<double>> Field_Map::find_fits(int position[3]){
    std::vector<std::vector<double>> abc_s;
    // Check dimensionality
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
	if (third_D_abc[0]>95.9728 && third_D_abc[0]<95.973){
	    double val_1 = points[6].get_value();
	    double val_2 = points[7].get_value();
	    double val_3 = points[8].get_value();
	    double diff_12 = val_1-val_2;
	    double diff_23 = val_2-val_3;
	    double diff_13 = val_1-val_3;
	    std::cout<<"FOUND:"<<" "<<position[0]<<","<<position[1]<<","<<position[2]<<"("<<shifted_three_d_central_pos<<")"<<std::endl;
	    std::cout<<"FOUND:"<<" "<<val_1<<std::endl;
	    std::cout<<"FOUND:"<<" "<<val_2<<std::endl;
	    std::cout<<"FOUND:"<<" "<<val_3<<std::endl;
	    std::cout<<"Diff 12:"<<" "<<diff_12<<std::endl;
	    std::cout<<"Diff 23:"<<" "<<diff_23<<std::endl;
	    std::cout<<"Diff 13:"<<" "<<diff_13<<std::endl;
	}
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
	// out of low bound: state edge but with warning
        closest_positions[0]=0;
	closest_positions[1]=0;
        _num_times_outside_domain[dimension-1][0]+=1;
    } else if (value>value_holder[0][2]){
	// out of high bound: state edge but with warning
        closest_positions[0]=pos_holder[0][2];
	closest_positions[1]=pos_holder[0][2];
	_num_times_outside_domain[dimension-1][1]+=1;
    } else{
    int test_point = 1;
    int num_secs = pos_holder.size();
    //std::cout<<"Num Secs: "<<num_secs<<std::endl;
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

threevector Field_Map::find_approx_value(threevector point, bool is_print){
    // Note to self: this entire method needs correcting to ensure that there is no excessive use of if statements and to ensure the distance is logical in spherical and cylindrical contexts.

    //std::cout<<"Started looking for value."<<std::endl;
    std::vector<int> one_d_close_pos(2);
    std::vector<int> two_d_close_pos(2);
    std::vector<int> three_d_close_pos(2);
    // number of points for 1 dimension = 2; 2d=4;3d=8
    const int num_near_points= std::pow(2, _dimension_val);
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
	    //std::cout<<"1D low: "<<one_d_close_pos[0]<<", 1D High: "<<one_d_close_pos[1]<<std::endl;
	    //std::cout<<"2D low: "<<two_d_close_pos[0]<<", 2D High: "<<two_d_close_pos[1]<<std::endl;
	    //std::cout<<"3D low: "<<three_d_close_pos[0]<<", 3D High: "<<three_d_close_pos[1]<<std::endl;
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

    //for averaging, linear but the closer the better
    double summed_weight = 0;
    threevector E_field_values[num_near_points];
    double weights[num_near_points];
    threevector summed_weighted_values = threevector(0,0,0);
    bool is_exact_pos_found = false;
    threevector exact_E_field_val;
    for (int i=0;i<num_near_points;i++){
        if (_is_E_Field_defined_truth || !_CALCULATE_NEARBY_E_FIELD_ON_CURVES) {
	    //E_field_values[i] = near_points[i].get_electric_field_at_point();
	    //Note to self: below not working
            E_field_values[i] = near_points[i].get_nearby_E_field(point, _dimension_val, _is_cartesian, _is_spherical, _is_cylindrical, _is_E_Field_defined_truth);
        } else {
            E_field_values[i] = near_points[i].get_nearby_E_field(point, _dimension_val, _is_cartesian, _is_spherical, _is_cylindrical, _is_E_Field_defined_truth);
        }
	if (is_print){
	    if (_is_cartesian) {
		std::cout<<i<<". Pos:    (x,y,z) "<< near_points[i].get_position_x()<<", "<< near_points[i].get_position_y()<<", "<<near_points[i].get_position_z()<<std::endl;
    	        std::cout<<i<<". E field:(x,y,z) "<<E_field_values[i].getx()<<", "<< E_field_values[i].gety()<<", "<<E_field_values[i].getz()<<std::endl;
	        std::cout<<"--"<<std::endl;
	    }
	    else if (_is_spherical) {
		std::cout<<i<<". Pos:    (r,theta,phi) "<< near_points[i].get_position_r()<<", "<< near_points[i].get_position_theta()<<", "<<near_points[i].get_position_phi()<<std::endl;
    	        std::cout<<i<<". E field:(r,theta,phi) "<<E_field_values[i].mag3()<<", "<< E_field_values[i].gettheta()<<", "<<E_field_values[i].getphi()<<std::endl;
	        std::cout<<"--"<<std::endl;
	    } else if (_is_cylindrical) {
		std::cout<<i<<". Pos:    (rho,z,phi) "<< near_points[i].get_position_rho()<<", "<< near_points[i].get_position_z()<<", "<<near_points[i].get_position_phi()<<std::endl;
    	        std::cout<<i<<". E field:(rho,z,phi) "<<E_field_values[i].getrho()<<", "<< E_field_values[i].getz()<<", "<<E_field_values[i].getphi()<<std::endl;
	        std::cout<<"--"<<std::endl;
	    }
	}
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
    if (is_exact_pos_found){
        return exact_E_field_val;
    }
    else {
        return threevector(summed_weighted_values.getx()/summed_weight, summed_weighted_values.gety()/summed_weight, summed_weighted_values.getz()/summed_weight);
    }
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
    const std::string symmetric_warning_message = "Warning! The user-defined" + coordinate_name + "range is not complete. Proceeding using extrapolation, but take caution.";
    if (dimension_to_check == 2){
        min = _two_d_value_holder[0][0];
        max = _two_d_value_holder[0][2];
        if (min > delta){
            const std::string min_string = std::to_str(min);
		    std::cout<<symmetric_warning_message<<"\n\tDefined lower limit: "<<min_string<<std::endl;
	    }
	    if (max < PI-delta){
            const std::string max_string = std::to_str(max);
		    std::cout<<symmetric_warning_message<<"\n\tDefined upper limit: "<<max_string<<std::endl;
	    }
    } else if (dimension_to_check == 3){
        min = _three_d_value_holder[0][0];
        max = _three_d_value_holder[0][2];
        if (min > delta){
            const std::string min_string = std::to_str(min);
		    std::cout<<symmetric_warning_message<<"\n\tDefined lower limit: "<<min_string<<std::endl;
	    }
	    if (max < 2*PI-delta){
            const std::string max_string = std::to_str(max);
		    std::cout<<symmetric_warning_message<<"\n\tDefined upper limit: "<<max_string<<std::endl;
	    }
    }
}

std::vector<double> Field_Map::find_complicit_domain_limits(bool is_spherical_injection, bool is_cylindrical_injection, double z_min_dimpl_pars, double z_max_dimpl_pars, double max_radius_dimpl_pars){
    double map_z_min;
    double map_z_max;
    double new_z_min;
    double new_z_max;
    const std::string z_limit_warning_message = "Warning! The user-defined custom potential map does not cover the full range of the inputted DiMPl pars. The z limit of the DiMPl pars shall be altered to the limit of the custom potential map.";

    double map_rho_max;
    double new_rho_max;
    const std::string rho_limit_warning_message = "Warning! The user-defined custom potential map does not cover the full range of the inputted DiMPl pars. The radial limit of the DiMPl pars shall be altered to the limit of the custom potential map.";

    if (_is_cartesian){
        map_z_min = _three_d_value_holder[0][0];
        map_z_max = _three_d_value_holder[0][2];
        const double limiting_x = std::min(std::abs(_one_d_value_holder[0][0]), _one_d_value_holder[0][2]);
        const double limiting_y = std::min(std::abs(_two_d_value_holder[0][0]), _two_d_value_holder[0][2]);
        map_rho_max = std::sqrt(limiting_x*limiting_x+limiting_y*limiting_y);
    } else if (_is_spherical) {
        const double map_r = _one_d_value_holder[0][2];
        if (is_spherical_injection){
            map_z_max = map_r;
            map_z_min = -map_r;
            map_rho_max = map_r;
        }
        if (is_cylindrical_injection){
            if (map_r < max_radius_dimpl_pars){
                // compress the radius of injection.
                const double limiting_z_dimpl_pars = std::min(std::abs(z_min_dimpl_pars), z_max_dimpl_pars);
                if (map_r < limiting_z_dimpl_pars){
                    // map is fully within dimpl pars, compress both rho and z to get as much coverage as possible
                    map_z_max = map_r/std::sqrt(3);
                    map_z_min = -map_z_max;
                    map_rho_max = map_r*std::sqrt(2/3);
                } else {
                    // map has coverage over some z limit so keep as much z as possible and compress rho accordingly
                    map_z_max = limiting_z_dimpl_pars;
                    map_z_min = -map_z_max;
                    map_rho_max = std::sqrt(map_r*map_r - limiting_z_dimpl_pars*limiting_z_dimpl_pars);
                }
            } else {
                // take the radius of injection to be incompressible.
                map_z_max = std::sqrt(map_r*map_r - max_radius_dimpl_pars*max_radius_dimpl_pars);
                map_z_min = -map_z_max;
                map_rho_max = map_r;
            }
        }
    } else if (_is_cylindrical){
        map_z_min = _two_d_value_holder[0][0];
        map_z_max = _two_d_value_holder[0][2];
        map_rho_max = _one_d_value_holder[0][2];
    }


    if (map_z_min > z_min_dimpl_pars){
        const std::string z_min_dimpl_pars_string = std::to_string(z_min_dimpl_pars);
        const std::string map_z_min_string = std::to_string(map_z_min);
        std::cout<<z_limit_warning_message<<"\n\t Old z lower bound (zmincoeff): "<<z_min_dimpl_pars<<"\n\t New z lower bound (zmincoeff): "<<map_z_min_string <<std::endl;
        new_z_min = map_z_min;
    } else {
        new_z_min = z_min_dimpl_pars; // keep old value
    }
    if (map_z_max < z_max_dimpl_pars){
        const std::string z_min_dimpl_pars_string = std::to_string(z_min_dimpl_pars);
        const std::string map_z_min_string = std::to_string(map_z_min);
        std::cout<<z_limit_warning_message<<"\n\t Old z Upper bound (zmaxcoeff): "<<z_min_dimpl_pars<<"\n\t New z upper bound (zmaxcoeff): "<<map_z_min_string <<std::endl;
        new_z_max = map_z_max;
    } else {
        new_z_max = z_max_dimpl_pars; // keep old value
    }

    if (map_rho_max < max_radius_dimpl_pars){
        const std::string max_radius_dimpl_pars_string = std::to_string(max_radius_dimpl_pars);
        const std::string map_rho_max_string = std::to_string(map_rho_max);
        std::cout<<rho_limit_warning_message<<"\n\t Old radial limit (impactpar): "<<max_radius_dimpl_pars<<"\n\t New radial limit (impactpar): "<<map_rho_max_string <<std::endl;
        new_rho_max = map_rho_max;
    } else {
        new_rho_max  = max_radius_dimpl_pars; // keep old value
    }

    std::vector<double> new_lim_vals = {new_z_min, new_z_max, new_rho_max};
    return new_lim_vals;
}

void Field_Map::check_dust_grain_map_coverage(double a1, double a2, double a3){
    double map_r_min;
    // Cartesian should be fine, so long as it is symmetric and covers dimpl domain which is checked elsewhere.

    if (_is_spherical || _is_cylindrical) {
        const double limiting_a = std::min(std::min(a1, a2), a3);
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
        if (dimension==1){coordinate_name = 'x';}
        else if (dimension==2){coordinate_name = 'y';}
        else if (dimension==3){coordinate_name = 'z';}
    } else if (_is_spherical){
        if (dimension==1){coordinate_name = 'r';}
        else if (dimension==2){coordinate_name = 'theta';}
        else if (dimension==3){coordinate_name = 'phi';}
    } else if (_is_cylindrical){
        if (dimension==1){coordinate_name = 'rho';}
        else if (dimension==2){coordinate_name = 'z';}
        else if (dimension==3){coordinate_name = 'phi';}
    }
    return coordinate_name;
}

void Field_Map::check_domain(bool is_spherical_injection, bool is_cylindrical_injection){
   // Check for symmetry and completeness.
   // Check for logical choise of coordinate system.
   // Check the z and radial injection limits set by dimpl are covered by the custom map
   // Check for dust limits
   const std::string conflicting_map_injection_warning_message = "Warning! The injection method coordinate system chosen is different to the custom potential map coordinate system.\n\tProceeding, but please be mindful.";

   if (_is_cartesian){
       if (is_spherical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tSpherical injection, cartesian potential map."<<std::endl;
       }
       if (is_cylindrical_injection){
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
       if (is_cylindrical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tCylindrical injection, spherical potential map."<<std::endl;
       }
       if (_dimension_val==2||_dimension_val==3){
           check_coordinate_domain_completeness(2);
           if (_dimension_val==3){
               check_coordinate_domain_completeness(3);
           }
       }
   } else if (_is_cylindrical){
       if (is_spherical_injection){
           std::cout<<conflicting_map_injection_warning_message<<"\n\tSpherical injection, cylindrical potential map."<<std::endl;
       }
       if (_dimension_val==2||_dimension_val==3){
           check_coordinate_domain_symmetry(2);
           if (_dimension_val==3){
               check_coordinate_domain_completeness(3);
           }
       }
   }
   check_dust_grain_map_coverage(a1, a2, a3);
}