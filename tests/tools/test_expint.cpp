
// Include general usage libraries
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../containers.hpp"
#include "../../src/config.hpp"
#include "../../src/special_math.hpp"
#include "../../src/tools.hpp"


bool launch_test(std::string file_path, std::function<cusfloat(cusfloat)> f_test)
{
    std::cout << "launching test..." << std::endl;
    // Read reference data
    DataRef data_ref;
    data_ref.read_single_channel(file_path);

    // Loop over reference data
    bool pass = true;
    cusfloat diff = 0.0;
    cusfloat fi = 0.0;
    std::cout << "Num.Points: " << data_ref.num_points << std::endl;
    for (int i=0; i<data_ref.num_points; i++)
    {
        // Calculate function value
        fi = f_test(data_ref.x[i]);

        // Check precision
        diff = std::abs(fi-data_ref.y[i]);
        if ((std::log10(diff)-std::log10(fi)) > EPS_PRECISION_ORDER)
        {
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

    std::string file_path_expinti(argv[1]);

    // Declare local variables
    bool pass = false;

    // Test exponential integral Ei(x)
    pass = launch_test(file_path_expinti, expint_i);
    if (!pass)
    {
        std::cerr << "test_expint/expint_i failed!" << std::endl;
        return 1;
    }

    return 0;
}