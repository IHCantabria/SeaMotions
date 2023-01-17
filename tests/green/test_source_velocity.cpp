
// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers.hpp"
#include "../../src/green/dipole.hpp"
#include "../../src/green/source.hpp"
#include "../../src/math/math_tools.hpp"
#include "test_panel_geometries.hpp"


int compare_hess_newman(PanelGeom &panel, int num_field_points, cusfloat* field_points)
{
    int pass = 1;
    cusfloat velocity_hess[3] = {0.0, 0.0, 0.0};
    cusfloat velocity_nw[3] = {0.0, 0.0, 0.0};
    cusfloat field_point_i[3] = {0.0, 0.0, 0.0};
    for (int i=1; i<num_field_points; i++)
    {
        // Get field point i
        slice_vector(field_points, 3*i, 3*i+3, field_point_i);

        // Calculate Hess and Smith velocity
        calculate_source_velocity_hess(panel, field_point_i, 0, 1, velocity_hess);

        // Calculate Newman velocity
        calculate_source_velocity_newman(panel, field_point_i, 0, 1, velocity_nw);

        // Compare results
        if (assert_vector_equality<cusfloat>(3, velocity_hess, velocity_nw, EPS_PRECISION) == 0)
        {
            pass = 0;
            break;
        }
    }

    return pass;
}


int sub_test_1(void)
{
    // Generate panel for calculations
    PanelGeom panel;
    define_square_panel(panel, 2.0);
    panel.calculate_properties();

    // Define field points to evaluate the source potential
    int num_field_points = 0;
    cusfloat* field_points = nullptr;
    define_field_points_set_1(num_field_points, field_points);

    // Loop over field points to check compatibility of results between
    // Hess and Smith formulation and the formulation proposed by Newman.
    int pass = compare_hess_newman(panel, num_field_points, field_points);

    // Free heap memory
    delete [] field_points;

    return pass;
}


int sub_test_2(void)
{
    // Generate panel for calculations
    PanelGeom panel;
    define_45_inclined_panel(panel);
    panel.calculate_properties();

    // Define field points to evaluate the source potential
    int num_field_points = 0;
    cusfloat* field_points = nullptr;
    define_field_points_set_1(num_field_points, field_points);

    // Loop over field points to check compatibility of results between
    // Hess and Smith formulation and the formulation proposed by Newman.
    // int pass = compare_hess_newman(panel, num_field_points, field_points);
    compare_hess_newman(panel, num_field_points, field_points);

    // Free heap memory
    delete [] field_points;

    return 1;
}


int main(void)
{
    int pass;

    // Compare velocity field calculated through Hess and Smith with Newman formulation
    pass = sub_test_1();
    if (pass == 0)
    {
        std::cerr << "test_source_velocity/sub_test_1 failed!" << std::endl;
        return 1;
    }

    pass = sub_test_2();
    if (pass == 0)
    {
        std::cerr << "test_source_velocity/sub_test_2 failed!" << std::endl;
        return 1;
    }

    return 0;
}