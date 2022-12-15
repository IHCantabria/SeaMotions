
// Include general usage libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Include local modules
#include "../containers.hpp"
#include "../../src/config.hpp"
#include "../../src/special_math.hpp"
#include "../../src/tools.hpp"


int main(int argc, char* argv[])
{
    // Check input arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path(argv[1]);

    // Read reference data for the test
    DataRef data_ref;
    data_ref.read_data(file_path);

    // Loop over the reference data values and check if the
    // output of the PSI function is correct according with the
    // machine precision defined during compilation
    cusfloat fi = 0.0;
    for (int i=0; i<data_ref.num_points; i++)
    {
        fi = psi_fun(static_cast<cusfloat>(data_ref.x[i]));
        if (std::abs(fi-data_ref.y[i]) > EPS_PRECISION)
        {
            std::cerr << "test_psi_function failed!" << std::endl;
        }
    }

    return 0;
}