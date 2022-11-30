
// Include local modules
#include "../../src/green/dipole.hpp"
#include "../../src/green/source.hpp"


void define_square_vertexes(PanelGeom &panel)
{
    // Define number of nodes for the panel
    panel.num_nodes = 4;

    // Define X position
    panel.x[0] = -1;
    panel.x[1] = -1;
    panel.x[2] = 1;
    panel.x[3] = 1;

    // Define Y position
    panel.y[0] = -1;
    panel.y[1] = 1;
    panel.y[2] = 1;
    panel.y[3] = -1;

    // Define Z position
    panel.z[0] = 0.0;
    panel.z[1] = 0.0;
    panel.z[2] = 0.0;
    panel.z[3] = 0.0;
}


void define_square_vertexes_2(PanelGeom &panel)
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


int main(void)
{
    // Define field point position
    cusfloat field_point[3] = {1.0, 1.0, std::sqrt(2.0)};
    // cusfloat field_point[3] = {1.0, 1.0, 1.0+std::sqrt(2)/2.0};

    // Define vertexes position
    PanelGeom panel;
    define_square_vertexes(panel);
    // define_square_vertexes_2(panel);
    panel.calculate_properties();

    // Calculate source velocity
    cusfloat velocity_nw[3] = {0.0, 0.0, 0.0};
    cusfloat velocity_hs[3] = {0.0, 0.0, 0.0};
    cusfloat phi = 0.0;
    calculate_source_velocity(panel, field_point, 0, velocity_nw);
    calculate_source_velocity_hess(panel, field_point, 0, velocity_hs);
    calculate_dipole_potential(panel, field_point, 0, phi);

    std::cout << "Velocity NW: "; print_vector(3, velocity_nw, 0);
    std::cout << "Velocity HS: "; print_vector(3, velocity_hs, 0);
    std::cout << "Phi: " << phi << std::endl;


    return 0;
}