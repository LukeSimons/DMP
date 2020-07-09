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
    std::cout<<"Opened constructor"<<std::endl;
    text_to_string_vector(field_file_name);
    std::cout<<"text to string vector complete"<<std::endl;
    find_header_details();
    std::cout<<"Header details found"<<std::endl;
    convert_data_line_to_multi_D_array();
    std::cout<<"Convereted data to multi d array"<<std::endl;
}

void Field_Map::convert_data_line_to_multi_D_array(){
    // Begin by finding the size
    _num_data_points_per_line = _dimension_val + 1;
    _num_data_lines = _line_vector.size() - _num_lines_in_header;
    
    // Create two D grid representation of data
    double two_D_array[_num_data_points_per_line][_num_data_lines];
    for (int i=0; i<_num_data_lines; i++){
	std::string line = _line_vector[i+_num_lines_in_header];
	std::string line_data[_num_data_points_per_line];
	for (int j=0; j<_num_data_points_per_line; j++){
	    int pos = line.find(_data_delim);
	    line_data[j] = line.substr(0, pos);
	    line.erase(0, pos + _data_delim.length());
	    two_D_array[j][i] = atof(line_data[j].c_str());
	}
    }
    // Convert grid representation to multi dimensional array of Field Points
    // Note to self: only working for three dimensional
    //
    // Begin by finding number of points in each dimension
    int one_d_count = 1;
    int two_d_count = 1;
    int three_d_count = 1;
    double memory_holder = two_D_array[1][0]; // start off with the first value
    // Note to self: this method only works if ordered
    for (int i=1; i< _num_data_lines; i++){
        if (two_D_array[1][i] < memory_holder){
	    break;
	}
	else {
	    memory_holder = two_D_array[1][i];
	    one_d_count+=1;
	}
    }
    memory_holder = two_D_array[2][0]; // start off with the first value
    for (int i=1; i<_num_data_lines/one_d_count; i++){
        if (two_D_array[2][i*one_d_count]<memory_holder){break;}
	else {
	    memory_holder = two_D_array[2][i*one_d_count];
	    two_d_count +=1;
	}
    }
    memory_holder = two_D_array[3][0]; // start off with the first value
    for (int i=1; i<_num_data_lines/one_d_count; i++){
        if (two_D_array[3][i*one_d_count*two_d_count]<memory_holder){break;}
	else {
	    memory_holder = two_D_array[3][i*one_d_count*two_d_count];
	    three_d_count +=1;
	}
    }


    // Populate multi dimensional array
    std::vector<std::vector<std::vector<Field_Point>>> Field_Map_Matrix(one_d_count, std::vector<std::vector<Field_Point>>(two_d_count, std::vector<Field_Point>(three_d_count)));
    for (int i=0; i<one_d_count; i++){
        for (int j=0; j<two_d_count; j++) {
	    for (int k=0; k<three_d_count; k++){
	        int two_d_array_pos = (k*one_d_count*two_d_count) + j * one_d_count + i; // find corresponding position in the two d grid
	        Field_Map_Matrix[i][j][k] = Field_Point(two_D_array[0][two_d_array_pos], threevector(two_D_array[1][two_d_array_pos], two_D_array[2][two_d_array_pos], two_D_array[3][two_d_array_pos])); // Note to self: need to include the type here
	    }
	}
    }
    _Field_Map_Final = Field_Map_Matrix;
}


void Field_Map::text_to_string_vector(std::string field_file_name)
{
    std::ifstream in(field_file_name);
    std::string single_line;
    std::vector<std::string> line_vector; // used to store contents of each line, vector used so that file can be arbitrary size
    
    if (in){
        while(getline(in, single_line))// go through each line in the file and write as a string
	{
	    _line_vector.push_back(single_line); // add each line as a new vector element
	}
    }
}

void Field_Map::find_header_details()
{
    _summary_line = split_after_delim(_summary_delim, _line_vector[0]);
    _dimension_val = std::stoi(split_after_delim(_dimension_delim, _line_vector[1])); // and convert string to integer
    _coord_type = split_after_delim(_coord_delim, _line_vector[2]);
    _ordered_truth = string_to_bool(split_after_delim(_ordered_delim, _line_vector[3])); // and convert string to boolean 
}

std::string Field_Map::split_after_delim(std::string delim, std::string str_to_split)
{
    return str_to_split.substr(str_to_split.find(delim) + delim.length(), str_to_split.size());
}

bool Field_Map::string_to_bool(std::string string){
    const std::string unknown_boolean_warning = "WARNING: ordered truth of custom field file is unknown, assuming ordered.";
	
    if (string.size() == 1){
        // the number case
	if (string.compare("0")){
	    return false;
	}
	else if (string.compare("1")){
	    return true;
	}
	else {
	    std::cout<< unknown_boolean_warning <<std::endl;
	    return true;
	}
    } else if (string.size() == 4){
        // the text boolean case
        std::transform(string.begin(), string.end(), string.begin(), [](unsigned char c) {return std::tolower(c); }); // convert string to lower case
	if (string.compare("true")){
	    return true;
	}
	else if (string.compare("false")){
	    return false;	
	}
        else {
	    std::cout << unknown_boolean_warning << std::endl;
	    return true;
	}
    }
    else {
        // unkown case
	std::cout << unknown_boolean_warning << std::endl;
	return true;
    }
}
