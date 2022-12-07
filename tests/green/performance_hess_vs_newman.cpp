
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers.hpp"
#include "../../src/green/source.hpp"
#include "test_panel_geometries.hpp"
#include "../../src/tools.hpp"


int main(void)
{
    // Define test parameters
    int num_reps = 1e7;

    // Create panel container
    PanelGeom panel;

    // Fill in with the unit square geometry
    define_square_panel(panel, 1.0);

    // Calcualte goemtric properties of the panel
    panel.calculate_properties();

    // Define field point
    cusfloat field_point[3] = {0.0, 0.0, 1.0};

    // Define velocity container
    cusfloat velocity[3] = {0.0, 0.0, 0.0};

    // Calculate time for Hess and Smith formulas
    double t0_hess = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        calculate_source_velocity_hess(panel, field_point, 0, 1, velocity);
    }
    double t1_hess = get_cpu_time();

    // Calculate time for Newman formulas
    double t0_newman = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        calculate_source_velocity_newman(panel, field_point, 0, 1, velocity);
    }
    double t1_newman = get_cpu_time();

    // Print results
    std::cout << "Elapsed time Hess: " << (t1_hess-t0_hess) << std::endl;
    std::cout << "Elapsed time Newman: " << (t1_newman-t0_newman) << std::endl;

    return 0;
}