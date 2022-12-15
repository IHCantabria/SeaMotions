
// Include general usage libraries
#include <fstream>
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


// Set relative precision for Bessel functions.
cusfloat EPS_BESSEL = 1e-6;


bool launch_test(cusfloat (*f)(cusfloat), std::string file_path, cusfloat precision, bool rel_flag)
{
    // Read reference data
    DataRef data_ref;
    data_ref.read_data(file_path);

    // Loop over reference data to check the solution of the 
    // function
    bool pass = true;
    cusfloat diff = 0.0;
    for (int i=0; i<data_ref.num_points; i++)
    {
        if (rel_flag)
        {
            diff = std::abs((f(data_ref.x[i]) - data_ref.y[i])/data_ref.y[i]);
        }
        else
        {
            diff = std::abs(f(data_ref.x[i]) - data_ref.y[i]);
        }

        if (diff > precision)
        {
            std::cout << std::setprecision(15) << "X: " << data_ref.x[i] << " - Fi(x): " << f(data_ref.x[i]);
            std::cout << std::setprecision(15) << " - Y_ref(x): " << data_ref.y[i] << " - Diff: " << diff << std::endl; 
            pass = false;
            break;
        }
    }

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 10))
    {
        return 1;
    }

    std::string file_path_besselj0(argv[1]);
    std::string file_path_besselj1(argv[2]);
    std::string file_path_bessely0(argv[3]);
    std::string file_path_bessely1(argv[4]);
    std::string file_path_bessels0(argv[5]);
    std::string file_path_bessels1(argv[6]);
    std::string file_path_besseli0(argv[7]);
    std::string file_path_besseli1(argv[8]);
    std::string file_path_besselk0(argv[9]);
    std::string file_path_besselk1(argv[10]);

    // Declare local variables
    bool pass = false;

    // Test first kind first order Bessel function
    pass = launch_test(besselj0, file_path_besselj0, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselj0 failed!" << std::endl;
        return 1;
    }

    // Test first kind second order Bessel function
    pass = launch_test(besselj1, file_path_besselj1, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselj1 failed!" << std::endl;
        return 1;
    }

    // Test second kind first order Bessel function
    pass = launch_test(bessely0, file_path_bessely0, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/bessely0 failed!" << std::endl;
        return 1;
    }

    // Test second kind second order Bessel function
    pass = launch_test(bessely1, file_path_bessely1, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/bessely1 failed!" << std::endl;
        return 1;
    }

    // Test first order struve function
    pass = launch_test(struve0, file_path_bessels0, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/struve0 failed!" << std::endl;
        return 1;
    }

    // Test second order struve function
    pass = launch_test(struve1, file_path_bessels1, EPS_BESSEL, false);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/struve1 failed!" << std::endl;
        return 1;
    }

    // Test mofified first kind first order Bessel function
    pass = launch_test(besseli0, file_path_besseli0, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besseli0 failed!" << std::endl;
        return 1;
    }

    // Test modified first kind second order Bessel function
    pass = launch_test(besseli1, file_path_besseli1, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besseli1 failed!" << std::endl;
        return 1;
    }

    // Test modified second kind first order Bessel function
    pass = launch_test(besselk0, file_path_besselk0, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselk0 failed!" << std::endl;
        return 1;
    }

    // Test modified second kind second order Bessel function
    pass = launch_test(besselk1, file_path_besselk1, EPS_BESSEL, true);
    if (!pass)
    {
        std::cerr << "test_bessel_functions/besselk1 failed!" << std::endl;
        return 1;
    }


    return 0;
}