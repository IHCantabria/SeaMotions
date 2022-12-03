
// Include general usage libraries
#include <fstream>
#include <iostream>
#include <string>

// Include local modules
#include "../../src/math_tools.hpp"
#include "test_panel_geometries.hpp"


void define_45_inclined_panel(PanelGeom &panel)
{
    // Define number of nodes for the panel
    panel.num_nodes = 4;

    // Define X position
    panel.x[0] = -1.0/std::sqrt(2);
    panel.x[1] = -1.0/std::sqrt(2);
    panel.x[2] = 1.0/std::sqrt(2);
    panel.x[3] = 1.0/std::sqrt(2);

    // Define Y position
    panel.y[0] = -1;
    panel.y[1] = 1;
    panel.y[2] = 1;
    panel.y[3] = -1;

    // Define Z position
    panel.z[0] = 2.0/std::sqrt(2);
    panel.z[1] = 2.0/std::sqrt(2);
    panel.z[2] = 0.0;
    panel.z[3] = 0.0;
}


void define_field_points_set_1(int &num_points, cusfloat* &data)
{
    // Define number of points
    num_points = 111;

    // Allocate heap memory for the array
    data = new cusfloat[3*num_points];

    // Declare z variable to use along the function
    cusfloat z;

    // Define points at z=0 m
    z = 0.0;
    int count_point = 0;
    data[3*count_point] = 0.0;
    data[3*count_point+1] = 0.0;
    data[3*count_point+2] = z;
    count_point++;

    fill_circle_points(data, 0.5, z, count_point);
    fill_circle_points(data, 1.0, z, count_point);
    fill_circle_points(data, 1.5, z, count_point);

    // Define points at z=1.0 m
    z = 1.0;
    data[3*count_point] = 0.0;
    data[3*count_point+1] = 0.0;
    data[3*count_point+2] = z;
    count_point++;
    
    fill_circle_points(data, 0.5, z, count_point);
    fill_circle_points(data, 1.0, z, count_point);
    fill_circle_points(data, 1.5, z, count_point);

    // Define points at z=-1.0 m
    z = -1.0;
    data[3*count_point] = 0.0;
    data[3*count_point+1] = 0.0;
    data[3*count_point+2] = z;
    count_point++;
    
    fill_circle_points(data, 0.5, z, count_point);
    fill_circle_points(data, 1.0, z, count_point);
    fill_circle_points(data, 1.5, z, count_point);

}


void define_square_panel(PanelGeom &panel, cusfloat scale)
{
    // Define number of nodes for the panel
    panel.num_nodes = 4;

    // Define X position
    panel.x[0] = -0.5*scale;
    panel.x[1] = -0.5*scale;
    panel.x[2] = 0.5*scale;
    panel.x[3] = 0.5*scale;

    // Define Y position
    panel.y[0] = -0.5*scale;
    panel.y[1] = 0.5*scale;
    panel.y[2] = 0.5*scale;
    panel.y[3] = -0.5*scale;

    // Define Z position
    panel.z[0] = 0.0;
    panel.z[1] = 0.0;
    panel.z[2] = 0.0;
    panel.z[3] = 0.0;
}


void fill_circle_points(cusfloat* data, cusfloat r, cusfloat z, int &count_point)
{
    for (int i=0; i<12; i++)
    {
        data[3*count_point] = r*std::cos(i*PI/12);
        data[3*count_point+1] = r*std::sin(i*PI/12);
        data[3*count_point+2] = z;
        
        count_point++;
    }

}