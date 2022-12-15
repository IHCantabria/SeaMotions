
// Include general usage libraries
#include <functional>
#include <iomanip>
#include <iostream>

// Include local modules
#include "../containers.hpp"
#include "../../src/config.hpp"
#include "../../src/special_math.hpp"
#include "../../src/tools.hpp"

#ifdef SIMPLE_PREC
cusfloat EPS_LEGENDRE = 1e-6;
#else
cusfloat EPS_LEGENDRE = 1e-12;
#endif


bool launch_test(DataRef &data_ref, int order)
{
    // Loop over reference data to check the target
    // function
    cusfloat fi = 0;
    cusfloat diff = 0;
    bool pass = true;
    for (int i=0; i<data_ref.num_points; i++)
    {
        fi = legendre_poly_raw(order, data_ref.x[i]);
        diff = data_ref.data[i][order] - fi;
        if (std::abs(diff)>EPS_LEGENDRE)
        {
            std::cerr << std::setprecision(15) << "X: " << data_ref.x[i] << " - Yf: " << fi;
            std::cerr << std::setprecision(15) << " - Yref:" << data_ref.data[i][order];
            std::cerr << " - Diff: " << diff << std::endl;
            pass = false;
            break;
        }
    }

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path_legendre(argv[1]);

    // Define local variables
    bool pass = false;

    // Load reference data
    DataRef data_ref;
    data_ref.read_multiple_channel(file_path_legendre);

    // Launch test for zeroth order Legendre polynomial
    for (int i=0; i<data_ref.num_channels; i++)
    {
        pass = launch_test(data_ref, i);
        if (!pass)
        {
            std::cerr << "test_legendre_poly/legendre_poly Order 0 failed!" << std::endl;
            return 1;
        }

    }

    return 0;
}