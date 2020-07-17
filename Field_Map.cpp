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
    add_quadratic_fit_details_to_all_points();
    duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
    std::cout<<"Fitted all points with a quadratic fit, time: "<<duration<<std::endl;
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
        std::coout<<"WARNING: coordinate system is not recognised, assumed Cartesian."<<std::endl;
        _is_cartesian = true;
    }
}

std::string Field_Map::split_after_delim(std::string delim, std::string str_to_split)
{
    return str_to_split.substr(str_to_split.find(delim) + delim.length(), str_to_split.size());
}

bool Field_Map::string_to_bool(std::string string){
    const std::string unknown_boolean_warning = "WARNING: ordered truth of custom field file is unknown, assuming ordered.";
    bool string_as_boolean = true; // initially assume ordered
    if (string.size() == 1){
        // the number case
	    if (string.compare("0")){
	        string_as_boolean = false;
	    }
	    else if (string.compare("1")){
	        string_as_boolean = true;
	    }
	    else {
	        std::cout<< unknown_boolean_warning <<std::endl;
	        string_as_boolean = true;
	    }
    } else if (string.size() == 4){
        // the text boolean case
        std::transform(string.begin(), string.end(), string.begin(), [](unsigned char c) {return std::tolower(c); }); // convert string to lower case
	    if (string.compare("true")){string_as_boolean = true;}
	    else if (string.compare("false")){string_as_boolean = false;}
        else {
	        std::cout << unknown_boolean_warning << std::endl;
	        string_as_boolean = true;
	    }
    }
    else {
        // unkown case
	    std::cout << unknown_boolean_warning << std::endl;
	    string_as_boolean = true;
    }
    return string_as_boolean;
}


void Field_Map::convert_data_line_to_multi_D_array(){
    // Begin by finding the size
    _num_data_points_per_line = _dimension_val + 1;
    _num_data_lines = _line_vector.size() - _num_lines_in_header;

    // Create two D grid representation of data
    double two_D_array[4][_num_data_lines];
    for (int i=0; i<_num_data_lines; i++){
	    std::string line = _line_vector[i+_num_lines_in_header];
	    std::string line_data[_num_data_points_per_line];
	    for (int j=0; j<_num_data_points_per_line; j++){
	        int pos = line.find(_data_delim);
	        line_data[j] = line.substr(0, pos);
	        line.erase(0, pos + _data_delim.length());
	        two_D_array[j][i] = atof(line_data[j].c_str());
	    }
        if (_dimension_val==1){
            two_D_array[2][i] = 0.0;
            two_D_array[3][i] = 0.0;
        }
        if (_dimension_val==2){
            two_D_array[3][i] = 0.0;
        }
    }
    // Convert grid representation to multi dimensional array of Field Points
    // Note to self: only working for three dimensional
    //
    // Begin by finding number of points in each dimension
    int one_d_count = 1;
    int two_d_count = 1;
    int three_d_count = 1;
    double memory_holder = two_D_array[1][0]; // start off with the first value in the first dimension
    if (!_ordered_truth){std::cout<<"WARNING: method only currently works for ordered files. Proceeding, but be mindful of potential errors."<<std::endl;}
    for (int i=1; i< _num_data_lines+1; i++){
        const double first_d_current_val = two_D_array[1][i];
        if ( first_d_current_val < memory_holder){
	        break;
	    }
	    else {
	        memory_holder = first_d_current_val;
	        one_d_count+=1;
	    }
    }
    if (_dimension_val==2||_dimension_vall==3){
        memory_holder = two_D_array[2][0]; // start off with the first value
        for (int i=1; i<_num_data_lines/one_d_count+1; i++){
            const int step = i*one_d_count;
            const double second_d_current_val = two_D_array[2][step];
            if (second_d_current_val<memory_holder){break;}
	        else {
	            memory_holder = second_d_current_val;
	            two_d_count +=1;
	        }
        }
        if (_dimension_val==3){
        memory_holder = two_D_array[3][0]; // start off with the first value
            for (int i=1; i<_num_data_lines/(one_d_count*two_d_count)+1; i++){
                const int step = i*one_d_count*two_d_count;
                const double third_d_current_val = two_D_array[3][step];
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
    _total_one_d_count = one_d_count;
    _total_two_d_count = two_d_count;
    _total_three_d_count = three_d_count;


    // Populate multi dimensional array
    std::vector<std::vector<std::vector<Field_Point>>> Field_Map_Matrix(one_d_count, std::vector<std::vector<Field_Point>>(two_d_count, std::vector<Field_Point>(three_d_count)));
    for (int i=0; i<one_d_count; i++){
        for (int j=0; j<two_d_count; j++) {
	        for (int k=0; k<three_d_count; k++){
	            int two_d_array_pos = (k*one_d_count*two_d_count) + j * one_d_count + i; // find corresponding position in the two d grid
                if (_is_cartesian){
                    Field_Map_Matrix[i][j][k] = Field_Point(two_D_array[0][two_d_array_pos], threevector(two_D_array[1][two_d_array_pos], two_D_array[2][two_d_array_pos], two_D_array[3][two_d_array_pos]));
                }
                else if (_is_spherical) {
                    Field_Map_Matrix[i][j][k] = Field_Point(two_D_array[0][two_d_array_pos], threevector(two_D_array[1][two_d_array_pos], two_D_array[2][two_d_array_pos], two_D_array[3][two_d_array_pos]), "s");
                }
                else if (_is_cylindrical) {
                    Field_Map_Matrix[i][j][k] = Field_Point(two_D_array[0][two_d_array_pos], threevector(two_D_array[1][two_d_array_pos], two_D_array[2][two_d_array_pos], two_D_array[3][two_d_array_pos]), "c");
                }
	        }
	    }
    }
    _Field_Map_Final = Field_Map_Matrix;
}

void Field_Map::partition_three_D_grid(){
    // Note to self: check how many dimensions are required
    // vector of vectors, in one direction the partition values, in the other direction additional granularity
    const std::vector<std::vector<int>> base_pos_holder= {{0}};
    const std::vector<std::vector<double>> base_value_holder= {{0.0}};
    if (_dimension_val==1){
        find_per_dimension(1);
        find_per_dimension(2);
        find_per_dimension(3);
    } else if (_dimension_val==2){
        find_per_dimension(1);
        find_per_dimension(2);
        _three_d_pos_holder = base_pos_holder;
        _three_d_value_holder = base_value_holder;
    } else if (_dimension_val==3){
        find_per_dimension(1);
        _two_d_pos_holder = base_pos_holder;
        _two_d_value_holder = base_value_holder;
        _three_d_pos_holder = base_pos_holder;
        _three_d_value_holder = base_value_holder;
    } else {
        std::cout<<"WARNING: dimension value unknown."
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
	    int num = pos_holder[i-1].size()+ 2*i;
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
    for (int k=0; k<_total_three_d_count; k++){
        for (int j=0; j<_total_two_d_count; j++){
	        for (int i=0; i<_total_one_d_count; i++){
		        int pos[3] = {i,j,k};
		        std::vector<std::vector<double>> fits = find_fits(pos);
		        _Field_Map_Final[i][j][k].add_quadratic_fit(fits);
		        _Field_Map_Final[i][j][k].find_electric_field();
	        }
	    }
    }
}

std::vector<std::vector<double>> Field_Map::find_fits(int position[3]){
    // If an edge situation, repeat the point
    bool is_first_D_on_low_edge = false;
    bool is_first_D_on_high_edge = false;
    bool is_second_D_on_low_edge = false;
    bool is_second_D_on_high_edge = false;
    bool is_third_D_on_low_edge = false;
    bool is_third_D_on_high_edge = false;
    Field_Point points[7];
    points[0] =  _Field_Map_Final[position[0]][position[1]][position[2]];
    if (position[0]==0){
        points[1]=points[0];
        is_first_D_on_low_edge=true;
    }
    else { points[1] =  _Field_Map_Final[position[0]-1][position[1]][position[2]];}

    if (position[0]==_total_one_d_count-1){
        points[2]=points[0];
        is_first_D_on_high_edge=true;
    }
    else { points[2] =  _Field_Map_Final[position[0]+1][position[1]][position[2]];}

    if (position[1]==0){
        points[3]=points[0];
        is_second_D_on_low_edge=true;
    }
    else { points[3] =  _Field_Map_Final[position[0]][position[1]-1][position[2]];}

    if (position[1]==_total_two_d_count-1){
        points[4]=points[0];
        is_second_D_on_high_edge=true;
    }
    else {points[4] =  _Field_Map_Final[position[0]][position[1]+1][position[2]];}

    if (position[2]==0){
        points[5]=points[0];
        is_third_D_on_low_edge=true;
    }
    else {points[5] =  _Field_Map_Final[position[0]][position[1]][position[2]-1];}

    if (position[2]==_total_three_d_count-1){
        points[6]=points[0];
        is_third_D_on_high_edge=true;
    }
    else {points[6] =  _Field_Map_Final[position[0]][position[1]][position[2]+1];}

    if (_is_cartesian){
        double first_D_p0[2] = {points[1].get_position_x(), points[1].get_value()};
        double first_D_p1[2] = {points[0].get_position_x(), points[0].get_value()};
        double first_D_p2[2] = {points[2].get_position_x(), points[2].get_value()};
        double second_D_p0[2] = {points[3].get_position_y(), points[3].get_value()};
        double second_D_p1[2] = {points[0].get_position_y(), points[0].get_value()};
        double second_D_p2[2] = {points[4].get_position_y(), points[4].get_value()};
        double third_D_p0[2] = {points[5].get_position_z(), points[5].get_value()};
        double third_D_p1[2] = {points[0].get_position_z(), points[0].get_value()};
        double third_D_p2[2] = {points[6].get_position_z(), points[6].get_value()};
    } else if(_is_spherical){
        double first_D_p0[2] = {points[1].get_position_r(), points[1].get_value()};
        double first_D_p1[2] = {points[0].get_position_r(), points[0].get_value()};
        double first_D_p2[2] = {points[2].get_position_r(), points[2].get_value()};
        double second_D_p0[2] = {points[3].get_position_theta(), points[3].get_value()};
        double second_D_p1[2] = {points[0].get_position_theta(), points[0].get_value()};
        double second_D_p2[2] = {points[4].get_position_theta(), points[4].get_value()};
        double third_D_p0[2] = {points[5].get_position_phi(), points[5].get_value()};
        double third_D_p1[2] = {points[0].get_position_phi(), points[0].get_value()};
        double third_D_p2[2] = {points[6].get_position_phi(), points[6].get_value()};
    } else if(_is_cylindrical){
        double first_D_p0[2] = {points[1].get_position_rho(), points[1].get_value()};
        double first_D_p1[2] = {points[0].get_position_rho(), points[0].get_value()};
        double first_D_p2[2] = {points[2].get_position_rho(), points[2].get_value()};
        double second_D_p0[2] = {points[3].get_position_z(), points[3].get_value()};
        double second_D_p1[2] = {points[0].get_position_z(), points[0].get_value()};
        double second_D_p2[2] = {points[4].get_position_z(), points[4].get_value()};
        double third_D_p0[2] = {points[5].get_position_phi(), points[5].get_value()};
        double third_D_p1[2] = {points[0].get_position_phi(), points[0].get_value()};
        double third_D_p2[2] = {points[6].get_position_phi(), points[6].get_value()};
    } else {std::cout<<"WARNING: coordinate type not found"<<std::endl;}

    std::vector<double> first_D_abc = calc_quad_fit(first_D_p0, first_D_p1, first_D_p2, is_first_D_on_low_edge, is_first_D_on_high_edge);
    std::vector<double> second_D_abc = calc_quad_fit(second_D_p0, second_D_p1, second_D_p2, is_second_D_on_low_edge, is_second_D_on_high_edge);
    std::vector<double> third_D_abc = calc_quad_fit(third_D_p0, third_D_p1, third_D_p2, is_third_D_on_low_edge, is_third_D_on_high_edge);
    std::vector<std::vector<double>> abc_s = {first_D_abc, second_D_abc, third_D_abc};
    return abc_s;
}

std::vector<double> Field_Map::calc_quad_fit(double p1[2], double p2[2], double p3[2], bool _is_repeated_low_point, bool _is_repeated_high_point){
    // solve in the form B = AC, C is unknown, A is matrix
    // Quadratic of the form y = ax^2+bx +c
    // C = {a,b,c}
    // B = {y1,y2,y3}
    std::vector<double> Y(3);
    if (_is_repeated_low_point&&_is_repeated_high_point){
	// Straight line between points if all match
        Y = {0,0,p1[1]};
    }
    else if(_is_repeated_low_point){
	    // use differing two for a linear fit between them
	    const double b = (p2[1]-p3[1])/(p2[0]-p3[0]);
	    const double c = p2[1]-b*p2[0];
	    Y = {0,b,c};
    }
    else if(_is_repeated_high_point){
	// use differing two for a linear fit between them
	const double b = (p1[1]-p2[1])/(p1[0]-p2[0]);
	const double c = p1[1]-b*p1[0];
	Y = {0,b,c};
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

std::vector<int> Field_Map::find_closest(double value, std::vector<std::vector<int>> pos_holder, std::vector<std::vector<double>> value_holder){

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

threevector Field_Map::find_approx_value(threevector point){
    // Note to self: this entire method needs correcting to ensure that there is no excessive use of fi statements and to ensure the distance is logical in spherical and cylindrical contexts.
    std::clock_t start;
    double duration;
    start = std::clock();
    std::cout<<"Started looking for value."<<std::endl;
    std::vector<int> one_d_close_pos(2);
    std::vector<int> two_d_close_pos(2);
    std::vector<int> three_d_close_pos(2);
    // number of points for 1 dimension = 2; 2d=4;3d=8
    const int num_near_points= std::pow(2, _dimension_val);
    Field_Point near_points[num_near_points];
    near_points[0] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][three_d_close_pos[0]];
    if(_dimension_val==1){
         if (_is_cartesian){
             one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder);
         } else if (_is_spherical){
             one_d_close_pos = find_closest(point.mag3(), _one_d_pos_holder, _one_d_value_holder);
         } else if (_is_cylindrical){
             one_d_close_pos = find_closest(point.getrho(), _one_d_pos_holder, _one_d_value_holder);
         }
         near_points[0] = _Field_Map_Final[one_d_close_pos[0]][0][0];
         near_points[1] = _Field_Map_Final[one_d_close_pos[1]][0][0];
    } else if (_dimension_val==2){
        if (_is_cartesian){
            one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder);
            two_d_close_pos = find_closest(point.gety(), _two_d_pos_holder, _two_d_value_holder);
        } else if (_is_spherical){
            one_d_close_pos = find_closest(point.mag3(), _one_d_pos_holder, _one_d_value_holder);
            two_d_close_pos = find_closest(point.gettheta(), _two_d_pos_holder, _two_d_value_holder);
        } else if (_is_cylindrical){
            one_d_close_pos = find_closest(point.getrho(), _one_d_pos_holder, _one_d_value_holder);
            two_d_close_pos = find_closest(point.gettheta(), _two_d_pos_holder, _two_d_value_holder);
        }
        near_points[0] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[0]][0];
        near_points[1] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[0]][0];
        near_points[2] = _Field_Map_Final[one_d_close_pos[0]][two_d_close_pos[1]][0];
        near_points[3] = _Field_Map_Final[one_d_close_pos[1]][two_d_close_pos[1]][0];
    } else if (_dimension_val==3){
        if (_is_cartesian){
            one_d_close_pos = find_closest(point.getx(), _one_d_pos_holder, _one_d_value_holder);
            two_d_close_pos = find_closest(point.gety(), _two_d_pos_holder, _two_d_value_holder);
            three_d_close_pos = find_closest(point.getz(), _three_d_pos_holder, _three_d_value_holder);
        } else if (_is_spherical){
            one_d_close_pos = find_closest(point.mag3(), _one_d_pos_holder, _one_d_value_holder);
            two_d_close_pos = find_closest(point.gettheta(), _two_d_pos_holder, _two_d_value_holder);
            three_d_close_pos = find_closest(point.getphi(), _three_d_pos_holder, _three_d_value_holder);
        } else if (_is_cylindrical){
            one_d_close_pos = find_closest(point.getrho(), _one_d_pos_holder, _one_d_value_holder);
            two_d_close_pos = find_closest(point.gettheta(), _two_d_pos_holder, _two_d_value_holder);
            three_d_close_pos = find_closest(point.getphi(), _three_d_pos_holder, _three_d_value_holder);
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
        E_field_values[i] = near_points[i].find_nearby_E_field(point);
	//std::cout<<"Pos:(i,j,k)"<< near_points[i].get_position_x()<<", "<< near_points[i].get_position_y()<<", "<<near_points[i].get_position_z()<<std::endl;
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
                distance_squared = std::pow(near_points[i].get_position_r()-point.mag3(),2) + std::pow(near_points[i].get_position_theta()-point.theta(),2);
            } else if (_is_cylindrical){
                distance_squared = std::pow(near_points[i].get_position_rho()-point.getrho(),2) + std::pow(near_points[i].get_position_z()-point.getz(),2);
            }
        } else if (_dimension_val==3){
            if (_is_cartesian){
                distance_squared = std::pow(near_points[i].get_position_x()-point.getx(),2)+std::pow(near_points[i].get_position_y()-point.gety(),2)+std::pow(near_points[i].get_position_z()-point.getz(),2);
            } else if (_is_spherical) {
                distance_squared = std::pow(near_points[i].get_position_r()-point.mag3(),2) + std::pow(near_points[i].get_position_theta()-point.theta(),2) + std::pow(near_points[i].get_position_phi()-point.phi(),2);
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
        duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
        std::cout<<"Finished looking for value, time: "<< duration<<std::endl;
        return exact_E_field_val;
    }
    else {
        duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
        std::cout<<"Finished looking for value, time: "<< duration<<std::endl;
        return threevector(summed_weighted_values.getx()/summed_weight, summed_weighted_values.gety()/summed_weight, summed_weighted_values.getz()/summed_weight);
    }
    duration = (std::clock() - start)/(double) CLOCKS_PER_SEC;
    std::cout<<"Finished looking for value, time: "<< duration<<std::endl;
}
