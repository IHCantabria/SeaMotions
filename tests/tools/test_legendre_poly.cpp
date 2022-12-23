
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


bool launch_channel_test(DataRef &data_ref, std::function<cusfloat(int, cusfloat)> f_test, int order)
{
    // Loop over reference data to check the target
    // function
    cusfloat fi = 0;
    cusfloat diff = 0;
    cusfloat diff_order = 0.0;
    bool pass = true;
    for (int i=0; i<data_ref.num_points; i++)
    {
        // Calculate current function value and its difference
        fi = f_test(order, data_ref.x[i]);
        diff = (data_ref.data[i][order] - fi);

        // Check for zero values in order to not have inifite
        // values when using the log10 check
        if (std::abs(diff)<1e-16)
        {
            diff = 1e-16;
        }

        if (std::abs(fi)<1e-16)
        {
            fi = 1e-16;
        }

        // Check precision
        diff_order = std::log10(std::abs(fi+1e-16))-std::log10(std::abs(diff));
        if (diff_order<std::log10(EPS_LEGENDRE))
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


bool launch_test(std::string file_path, std::function<cusfloat(int, cusfloat)> f_test,
    std::string test_type)
{
    // Load reference data
    DataRef data_ref;
    data_ref.read_multiple_channel(file_path);

    // Launch test for zeroth order Legendre polynomial
    bool pass = false;
    for (int i=0; i<data_ref.num_channels; i++)
    {
        pass = launch_channel_test(data_ref, f_test, i);
        if (!pass)
        {
            std::cerr << "test_legendre_poly/" << test_type;
            std::cerr << " Order:" << i << " failed!" << std::endl;
            break;
        }

    }

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 2))
    {
        return 1;
    }

    std::string file_path_legendre(argv[1]);
    std::string file_path_legendre_der(argv[2]);

    // Define local variables
    bool pass = false;

    // Test legendre polynomial function
    pass = launch_test(file_path_legendre, legendre_poly_raw, "legendre_poly");
    if (!pass)
    {
        return 1;
    }

    // Test legendre derivative polynomial function
    pass = launch_test(file_path_legendre_der, legendre_poly_der_raw, "legendre_poly_der");
    if (!pass)
    {
        return 1;
    }

    return 0;
}