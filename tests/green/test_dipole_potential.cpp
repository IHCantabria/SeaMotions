
// Include external libraries
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/containers.hpp"
#include "../../src/green/dipole.hpp"
#include "../../src/math/math_interface.hpp"
#include "../../src/math/math_tools.hpp"
#include "test_panel_geometries.hpp"
#include "../../src/tools.hpp"


int sub_test_1(void)
{
    // Define field point position
    cusfloat field_point[3] = {0.0, 0.0, std::sqrt(2.0)};

    // Define vertexes position
    PanelGeom panel;
    define_square_panel(panel, 2.0);
    panel.calculate_properties();

    // Calculate dipole potential
    cusfloat phi = 0.0;
    calculate_dipole_potential(panel, field_point, 0, phi);

    // Calculate analytical solution
    cusfloat phi_ana = 4*std::atan(1/2.0/std::sqrt(2.0));

    int pass = 1;
    if (std::abs(phi-phi_ana) > EPS_PRECISION)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_2(void)
{
    // Create unit square panel
    PanelGeom panel;
    define_square_panel(panel, 2.0);
    panel.calculate_properties();

    // Define field points to evaluate the source potential
    int num_field_points = 0;
    cusfloat* field_points = nullptr;
    define_field_points_set_1(num_field_points, field_points);

    // Loop over the field points and check equivalence in between Hess and Smith formulation
    // and Newmann method.
    int pass = 1;
    cusfloat phi_hess = 0.0;
    cusfloat phi_nm = 0.0;
    // for (int i=0; i<num_field_points; i++)
    // {
        // Calculate Hess and Smith potential
        // calculate_dipole_potential()

        // Calculate Newman potential

        // Compare
    // }

    // Free heap memory
    delete [] field_points;

    return pass;
}


int main(void)
{
    int pass;

    // Calculate dipole potential for a field point over the unit square panel.
    pass = sub_test_1();
    if (pass == 0)
    {
        std::cerr << "test_dipole_potential/sub_test_1 failed!" << std::endl;
        return 1;
    }

    // Calculate compatibility of source potential of Hess and Smith and Newman formulations
    pass = sub_test_2();
    if (pass == 0)
    {
        std::cerr << "test_dipole_potential/sub_test_2 failed!" << std::endl;
    }

    return 0;
}