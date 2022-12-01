
// Include local modules
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