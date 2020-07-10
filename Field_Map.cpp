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
    partition_three_D_grid();
    std::cout<<"Partitioned three D grid"<<std::endl;
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

    _total_one_d_count = one_d_count;
    _total_two_d_count = two_d_count;
    _total_three_d_count = three_d_count;


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

void Field_Map::partition_three_D_grid(){
    // Note to self: check how many dimensions are required
    // vector of vectors, in one direction the partition values, in the other direction additional granularity
    // Check what precision is required

    // Note to self: does this work on non-integer half values?
    //     aim to address this is by always rounding down and having an additional partition at the end
    // First Dimension
    std::vector<std::vector<double>> value_holder;
    std::vector<std::vector<int>> pos_holder;
    int num_secs = log2(_total_one_d_count);// find the log_2 value
    std::vector<double> line_one_val_holder(3); //line 1
    std::vector<int> line_one_pos_holder = partition(0, _total_one_d_count-1); // line 1
    for (int i=0; i<3; i++){
	line_one_val_holder[i] = _Field_Map_Final[line_one_pos_holder[i]][0][0].get_position_x();
    }
    pos_holder.push_back(line_one_pos_holder);
    value_holder.push_back(line_one_val_holder);
    for (int i=1; i<num_secs; i++){
	int num = pos_holder[i-1].size()+ 2;
	if (num>_total_one_d_count/3){
	    // need to go every single 1, this is the last case
	    num = _total_one_d_count;
	    std::cout<<"num:"<<num<<std::endl;
	}
	std::vector<double> this_iteration_val_holder(num);
	std::vector<int> this_iteration_pos_holder(num);
	if (num != _total_one_d_count){
	    this_iteration_pos_holder[0] = 0; // first index
	    for (int j=0; j<num-3; j++){
	        std::vector<int> next_three = partition(pos_holder[i-1][j], pos_holder[i-1][j+1]);
	        std::copy(next_three.begin(), next_three.end(), this_iteration_pos_holder.begin()+j*2);
	    }
	} else {
	    for (int j=0;j<_total_one_d_count;j++){
	        this_iteration_pos_holder[j]=j;
	    }
	}
	for (int j=0; j<num; j++){
	    this_iteration_val_holder[j] = _Field_Map_Final[this_iteration_pos_holder[j]][0][0].get_position_x();
	}
        value_holder.push_back(this_iteration_val_holder);
	pos_holder.push_back(this_iteration_pos_holder);
    }
    for (int i=0; i<num_secs; i++){
        std::cout<<"Row: ";
	int num_in_row = pos_holder[i].size();
	for (int j=0; j<num_in_row; j++){
	    std::cout << pos_holder[i][j]<<" ("<<value_holder[i][j]<<")" <<"; ";
	}
	std::cout<<std::endl;
    }
    find_closest(0.8, pos_holder, value_holder);
    find_closest(0.1, pos_holder, value_holder);
    find_closest(5, pos_holder, value_holder);
    find_closest(11.3, pos_holder, value_holder);
}

void Field_Map::find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder){
    std::cout<<"find closest called, looking for: "<<value<<std::endl;
    int test_point = 1;
    int num_secs = pos_holder.size();
    for (int i=0; i<num_secs-2; i++){
	if (value > value_holder[i][test_point]){
	    // occurs after midpoint, therefore check after here in the next line
            test_point+=2;
	}
    }
    // Penultimate line:
    //     check partition above or below depending on 
    std::vector<int> closest_positions(2);
    // Note to self: this can be improved by rahter than just finding the closest position, find the difference between the two closest points and evaluate accordingly
    std::cout<<"Value: "<<value<<"; test value: "<<value_holder[num_secs-2][test_point]<<std::endl;
    if (value < value_holder[num_secs-2][test_point]){
        test_point -= 1;
    }
        // check diff to this value, diff to next value etc. until buffer
        int num_to_check = 1+pos_holder[num_secs-2][test_point+1] - pos_holder[num_secs-2][test_point];
	std::cout<<"Number needed to check: "<<num_to_check<<std::endl;
	int index_of_test_point = pos_holder[num_secs-2][test_point];
	std::cout<<"Index of test point: "<<index_of_test_point<<std::endl;
        int closest_pos = index_of_test_point;
	double closest_diff = value-value_holder[num_secs-2][test_point];
	std::cout<<"Initial closest diff: "<< closest_diff<<std::endl;
	int second_closest_pos = index_of_test_point+1;
	for (int i=1; i<num_to_check;i++){
	    std::cout<<"Check location: "<<pos_holder[num_secs-1][index_of_test_point+i]<<std::endl;;
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
	std::cout<<"The closest position is: "<<closest_pos<<std::endl;
	std::cout<<"The second closest position is: "<<second_closest_pos<<std::endl;
        closest_positions[0] = closest_pos;
	closest_positions[1] = second_closest_pos;
    
    std::cout<<"Closest point: "<<closest_positions[0]<<" ("<<value_holder[num_secs-1][closest_positions[0]]<<")"<<std::endl;
    std::cout<<"Second Closest point: "<<closest_positions[1]<<" ("<<value_holder[num_secs-1][closest_positions[1]]<<")"<<std::endl;

}



std::vector<int> Field_Map::partition(int start, int end){
    int midway = start+(end-start)/2;
    return {start, midway, end};
}

double Field_Map::find_E_Val(threevector Position, double EField_Background){
    // Note to self: check to be in the same coordinate system for consistency
    return 0.0;    
}
