
// Include external libraries
#include <cstdlib>
#include <iostream>
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers.hpp"
#include "../../src/green/dipole.hpp"
#include "../../src/math_interface.hpp"
#include "../../src/math_tools.hpp"
#include "../../src/tools.hpp"


void define_square_vertexes(PanelGeom &panel)
{
    // Define number of nodes for the panel
    panel.num_nodes = 4;

    // Define X position
    panel.x[0] = -1;
    panel.x[1] = 1;
    panel.x[2] = 1;
    panel.x[3] = -1;

    // Define Y position
    panel.y[0] = -1;
    panel.y[1] = -1;
    panel.y[2] = 1;
    panel.y[3] = 1;

    // Define Z position
    panel.z[0] = 0.0;
    panel.z[1] = 0.0;
    panel.z[2] = 0.0;
    panel.z[3] = 0.0;
}


int main(void)
{
    // int num_points = 4;
    // double* x_pos = generate_empty_vector(num_points);
    // double* y_pos = generate_empty_vector(num_points);
    // double* z_pos = generate_empty_vector(num_points);
    // // Print points position
    // for (int i=0; i<10; i++)
    // {
    //     std::cout << x_pos[i] << " - " << y_pos[i] << " - " << z_pos[i] << std::endl;
    // }
    // // Free memory
    // mkl_free( x_pos );
    // mkl_free( y_pos );
    // mkl_free( z_pos );

    // Define field point position
    cusfloat field_point[3] = {0.0, 0.0, 1.0};

    // Define vertexes position
    PanelGeom panel;
    define_square_vertexes(panel);
    panel.calculate_properties();

    // double t0 = get_cpu_time();
    // for (int i=0; i<1e7; i++)
    // {
    //     panel.calculate_properties();
    // }
    // double t1 = get_cpu_time();
    // double elapsed = t1-t0;
    // std::cout << "Elapsed Time: " << elapsed << std::endl;

    // std::cout << "Precision: " << FLOATING_PRECISION << std::endl;
    // std::cout << "Floating type: " << typeid(cusfloat).name() << std::endl;
    // std::cout << "Floating size: " << sizeof(cusfloat) << std::endl;

    // std::cout << "PANEL NODES POSITION" << std::endl;
    // for (int i=0; i<4; i++)
    // {
    //     std::cout << "Point: " << i << std::endl;
    //     std::cout << "  X: " << panel.x[i] << std::endl;
    //     std::cout << "  Y: " << panel.y[i] << std::endl;
    //     std::cout << "  Z: " << panel.z[i] << std::endl;
    //     std::cout << std::endl;
    // }
    // std::cout << "Panel Center" << std::endl;
    // std::cout << "  X: " << panel.xm << std::endl;
    // std::cout << "  Y: " << panel.ym << std::endl;
    // std::cout << "  Z: " << panel.zm << std::endl;

    // Calculate dipole potential
    dipole_potential(panel, field_point);

    std::cout << "Program working... " << std::endl;

    return 0;
}