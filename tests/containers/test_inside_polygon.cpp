
// Include general usage libraries
#include <fstream>
#include <iostream>
#include <string>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/containers.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/tools.hpp"


void read_file_data(std::string file_path, PanelGeom &panel, int &num_field_points,
    cusfloat* &field_points, int* &is_inside)
{
    // Open file unit
    std::ifstream infile(file_path.c_str());

    // Read number of panel nodes
    infile >> panel.num_nodes;

    // Loop over number of nodes and extrac data from file. In this case the nodes 
    // value are pasted directly over the local nodes definition in order to not change 
    // from reference system and to have compatibility with the external test definition
    for (int i=0; i<panel.num_nodes; i++)
    {
        // Read from file
        infile >> panel.xl[i] >> panel.yl[i];

        // Clean-up memory just in case (not used)
        panel.zl[i] = 0.0;
    }

    // Read number of field points
    infile >> num_field_points;

    // Allocate heap memory for field points
    field_points = new cusfloat[3*num_field_points];
    is_inside = new int[num_field_points];

    // Read field points
    for (int i=0; i<num_field_points; i++)
    {
        // Read data from file
        infile >> field_points[3*i+0] >> field_points[3*i+1] >> is_inside[i];
        
        // Clean-up memory just in case
        field_points[3*i+2] = 0.0;
    }

    // Close file unit
    infile.close();
}


int test(std::string file_path)
{
    // Declare local variables
    int pass = 1;

    // Create panel container
    PanelGeom panel;

    // Create field points container
    int num_field_points = 0;
    cusfloat* field_points = nullptr; 

    // Create solution container
    int* is_inside = nullptr;

    // Read data from file
    read_file_data(file_path, panel, num_field_points, field_points, is_inside);

    // Loop over field points to check if they are inside the polygon or not
    int is_inside_i = 0;
    cusfloat field_point_aux[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<num_field_points; i++)
    {
        // Get field point i
        slice_vector(field_points, 3*i, 3*i+3, field_point_aux);

        // Check if it is inside
        is_inside_i = panel.is_inside(field_point_aux);

        // Compare calculation with the reference value
        if (is_inside_i != is_inside[i])
        {
            pass = 0;
            break;
        }
    }

    // Delete heap memory
    delete [] field_points;

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 3))
    {
        return 1;
    }
    
    std::string file_path_1(argv[1]);
    std::string file_path_2(argv[2]);
    std::string file_path_3(argv[3]);

    // Declare local variables
    int pass = 0;

    // Add sub-test 1
    pass = test(file_path_1);
    if (pass == 0)
    {
        std::cerr << "test_inside_polygon/sub_test_1 failed!" << std::endl;
        return 1;
    }

    // Add sub-test 2
    pass = test(file_path_2);
    if (pass == 0)
    {
        std::cerr << "test_inside_polygon/sub_test_2 failed!" << std::endl;
        return 1;
    }

    // Add sub-test 3
    pass = test(file_path_3);
    if (pass == 0)
    {
        std::cerr << "test_inside_polygon/sub_test_3 failed!" << std::endl;
        return 1;
    }

    return 0;
}