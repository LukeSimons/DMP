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
    int pos[3] = {3,2,1};
    std::vector<std::vector<double>> fits = find_fits(pos);
    std::cout<<"TEST=================="<<std::endl;
    Field_Point centre_point = _Field_Map_Final[pos[0]][pos[1]][pos[2]]; 
    std::cout<<"Centre Point Values (val, x,y,z): "<<centre_point.get_value()<<", "<<centre_point.get_position_x()<<", "<<centre_point.get_position_y()<<", "<< centre_point.get_position_z()<<std::endl; 
    std::cout<<"x: (a,b,c)"<<fits[0][0]<<", "<< fits[0][1]<<", "<<fits[0][2]<<std::endl;
    std::cout<<"y: (a,b,c)"<<fits[1][0]<<", "<< fits[1][1]<<", "<<fits[1][2]<<std::endl;
    std::cout<<"z: (a,b,c)"<<fits[2][0]<<", "<< fits[2][1]<<", "<<fits[2][2]<<std::endl;
    std::cout<<"====================="<<std::endl;
}

/** @brief Take raw data (in an ordered fashion) and transform it to a multi dimensional array.
 *  The multi dimensional array spans real space and holds a field value at each point.
 *  @param 
 *
 */
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

/** @brief Converts a text file to a vector containing contents of each line.
 *         Vector is of initially arbitrary size and each element is the line
 *         contents written as a string.
 *  @param field_file_name : the name of the text file
 *
 */
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
    find_per_dimension(1);
    find_per_dimension(2);
    find_per_dimension(3);

    std::vector<std::vector<int>> pos_holder = _two_d_pos_holder;
    std::vector<std::vector<double>> value_holder = _two_d_value_holder;
}

double Field_Map::find_approx_value(threevector point){
    std::vector<int> one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder);
    std::vector<int> two_d_close_pos = find_closest(point.gety(), _two_d_pos_holder, _two_d_value_holder);
    std::vector<int> three_d_close_pos = find_closest(point.getz(), _three_d_pos_holder, _three_d_value_holder);
    Field_Point near_points[8];
    near_points[0] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][three_d_close_pos[0]];
    near_points[1] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[0]][three_d_close_pos[0]];
    near_points[2] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[1]][three_d_close_pos[0]];
    near_points[3] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[1]][three_d_close_pos[0]];
    near_points[4] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][three_d_close_pos[1]];
    near_points[5] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[0]][three_d_close_pos[1]];
    near_points[6] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[1]][three_d_close_pos[1]];
    near_points[7] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[1]][three_d_close_pos[1]];
    //for averaging, linear but the closer the better
    double summed_weight = 0;
    double values[8];
    double weights[8];
    double summed_weighted_values = 0;
    bool is_exact_pos_found = false;
    double exact_val;
    for (int i=0;i<8;i++){
	values[i] = near_points[i].get_value();
	std::cout<<"Value: "<<values[i]<<std::endl;
	std::cout<<"Pos:(i,j,k)"<< near_points[i].get_position_x()<<", "<< near_points[i].get_position_y()<<", "<<near_points[i].get_position_z()<<std::endl;
        double distance_squared = std::pow(near_points[i].get_position_x()-point.getx(),2)+std::pow(near_points[i].get_position_y()-point.gety(),2)+std::pow(near_points[i].get_position_z()-point.getz(),2);
	if (distance_squared == 0){
	    is_exact_pos_found = true;
	    exact_val = values[i];
	}
	else {
	weights[i] = 1/distance_squared;
	summed_weight += weights[i];
	summed_weighted_values += values[i]*weights[i];
	}
    }
    if (is_exact_pos_found){return exact_val;}
    else {return summed_weighted_values/summed_weight;}
}

void Field_Map::find_per_dimension(int dimension){
    int count=1;//initialise with meaningless value
    std::string warning = "WARNING: invalid dimension provided";
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
    int num_secs = log2(count);// find the log_2 value
    std::vector<double> line_one_val_holder(3); //line 1
    std::vector<int> line_one_pos_holder = partition(0, count-1); // line 1
    if(dimension==1){
        for (int i=0; i<3; i++){
	    line_one_val_holder[i] = _Field_Map_Final[line_one_pos_holder[i]][0][0].get_position_x();
        }

    }else if (dimension==2){
        for (int i=0; i<3; i++){
	    line_one_val_holder[i] = _Field_Map_Final[0][line_one_pos_holder[i]][0].get_position_y();
        }

    } else if (dimension==3){
        for (int i=0; i<3; i++){
	    line_one_val_holder[i] = _Field_Map_Final[0][0][line_one_pos_holder[i]].get_position_z();
        }

    } else{
        std::cout<<warning<<std::endl;
    }
    pos_holder.push_back(line_one_pos_holder);

    value_holder.push_back(line_one_val_holder);

    for (int i=1; i<num_secs; i++){
	int num = pos_holder[i-1].size()+ 2*i;
	if (num>count/2){
	    // need to go every single 1, this is the last case
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
	    for (int j=0;j<count;j++){
	        this_iteration_pos_holder[j]=j;
	    }
	}

        if(dimension==1){
	    for (int j=0; j<num; j++){
	        this_iteration_val_holder[j] = _Field_Map_Final[this_iteration_pos_holder[j]][0][0].get_position_x();
	    }

        }else if (dimension==2){
	    for (int j=0; j<num; j++){
	        this_iteration_val_holder[j] = _Field_Map_Final[0][this_iteration_pos_holder[j]][0].get_position_y();
	    }

        } else if (dimension==3){
	    for (int j=0; j<num; j++){
	        this_iteration_val_holder[j] = _Field_Map_Final[0][0][this_iteration_pos_holder[j]].get_position_z();
	    }

        } else{
            std::cout<<warning<<std::endl;
        }

        value_holder.push_back(this_iteration_val_holder);
	pos_holder.push_back(this_iteration_pos_holder);
    }
    for (int i=0; i<num_secs; i++){
	int num_in_row = pos_holder[i].size();
	for (int j=0; j<num_in_row; j++){
	}
    }
    if(dimension==1){
        _one_d_pos_holder = pos_holder;
	_one_d_value_holder = value_holder;
    }else if (dimension==2){
        _two_d_pos_holder = pos_holder;
	_two_d_value_holder = value_holder;
    } else if (dimension==3){
        _three_d_pos_holder = pos_holder;
	_three_d_value_holder = value_holder;
    } else{
        std::cout<<warning<<std::endl;
    }
}

std::vector<int> Field_Map::find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder){
    /**
        Begin by looking at the values in each line.
        For the first line, if the true value is above the midpoint then on the line beneath the next midpoint should be considered.
        This is enacted by keeping track of a `test_point' which is shifted accordingly for each subsequent line.
        This is continued up until the line before the penultimate line (since the final line is spaced differently)
    */
    int test_point = 1;
    int num_secs = pos_holder.size();
    for (int i=0; i<num_secs-2; i++){
	      if (value > value_holder[i][test_point]){
            test_point+=std::pow(2,(i+1)); // for the next line, the test point is shifted
	      } else {
	          test_point+=test_point-1;
	      }
    }
    // Penultimate line:
    //     check partition above or below depending on
    std::vector<int> closest_positions(2);
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

    return closest_positions;
}



std::vector<int> Field_Map::partition(int start, int end){
    int midway = start+(end-start)/2;
    return {start, midway, end};
}

std::vector<std::vector<double>> Field_Map::find_fits(int position[3]){
    // If an edge situation, shift to the neighbouring one
    for (int i=0;i<3;i++){
    	if(position[i]==0){position[i]=1;}
    }
    if(position[0]==_total_one_d_count-1){position[0]=_total_one_d_count-2;}
    if(position[1]==_total_two_d_count-1){position[0]=_total_two_d_count-2;}
    if(position[2]==_total_three_d_count-1){position[0]=_total_three_d_count-2;}

    Field_Point points[7];
    points[0] =  _Field_Map_Final[position[0]][position[1]][position[2]];
    points[1] =  _Field_Map_Final[position[0]-1][position[1]][position[2]];
    points[2] =  _Field_Map_Final[position[0]+1][position[1]][position[2]];
    points[3] =  _Field_Map_Final[position[0]][position[1]-1][position[2]];
    points[4] =  _Field_Map_Final[position[0]][position[1]+1][position[2]];
    points[5] =  _Field_Map_Final[position[0]][position[1]][position[2]-1];
    points[6] =  _Field_Map_Final[position[0]][position[1]][position[2]+1];
    double first_D_p0[2] = {points[1].get_position_x(), points[1].get_value()};
    double first_D_p1[2] = {points[0].get_position_x(), points[0].get_value()};
    double first_D_p2[2] = {points[2].get_position_x(), points[2].get_value()};
    double second_D_p0[2] = {points[3].get_position_y(), points[3].get_value()};
    double second_D_p1[2] = {points[0].get_position_y(), points[0].get_value()};
    double second_D_p2[2] = {points[4].get_position_y(), points[4].get_value()};
    double third_D_p0[2] = {points[5].get_position_z(), points[5].get_value()};
    double third_D_p1[2] = {points[0].get_position_z(), points[0].get_value()};
    double third_D_p2[2] = {points[6].get_position_z(), points[6].get_value()};
    std::vector<double> first_D_abc = calc_quad_fit(first_D_p0, first_D_p1, first_D_p2);
    std::vector<double> second_D_abc = calc_quad_fit(second_D_p0, second_D_p1, second_D_p2);
    std::vector<double> third_D_abc = calc_quad_fit(third_D_p0, third_D_p1, third_D_p2);

    std::vector<std::vector<double>> abc_s = {first_D_abc, second_D_abc, third_D_abc};

    return abc_s;
}

std::vector<double> Field_Map::calc_quad_fit(double p1[2], double p2[2], double p3[2]){
    // solve in the form B = AC, C is unknown, A is matrix
    // Quadratic of the form y = ax^2+bx +c
    // C = {a,b,c}
    // B = {y1,y2,y3}
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
    std::vector<double>  Y = {y1,y2,y3};
    return Y;
}

double Field_Map::two_by_two_mat_det(double upper_left, double upper_right, double bottom_left, double bottom_right){
    return upper_left*bottom_right - upper_right*bottom_left;
}
