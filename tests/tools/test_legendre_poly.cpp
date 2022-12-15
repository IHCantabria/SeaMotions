
// Include general usage libraries
#include <functional>
#include <iostream>

// Include local modules
#include "../containers.hpp"
#include "../../src/config.hpp"
#include "../../src/special_math.hpp"
#include "../../src/tools.hpp"


bool launch_test(std::string file_path, int order)
{
    // Load reference data
    DataRef data_ref;
    data_ref.read_data(file_path);

    // Loop over reference data to check the target
    // function
    cusfloat fi = 0;
    cusfloat diff = 0;
    bool pass = true;
    for (int i=0; i<data_ref.num_points; i++)
    {
        fi = legendre_poly_raw(order, data_ref.x[i]);
        diff = data_ref.y[i] - fi;
        if (std::abs(fi-diff)>EPS_PRECISION)
        {
            std::cerr << "X: " << data_ref.x[i] << " - Yf: " << fi;
            std::cerr << " - Yref:" << data_ref.y[i] << std::endl;
            pass = false;
            break;
        }
    }

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 11))
    {
        return 1;
    }

    std::string file_path_legendre0(argv[1]);
    std::string file_path_legendre1(argv[2]);
    std::string file_path_legendre2(argv[3]);
    std::string file_path_legendre3(argv[4]);
    std::string file_path_legendre4(argv[5]);
    std::string file_path_legendre5(argv[6]);
    std::string file_path_legendre6(argv[7]);
    std::string file_path_legendre7(argv[8]);
    std::string file_path_legendre8(argv[9]);
    std::string file_path_legendre9(argv[10]);
    std::string file_path_legendre10(argv[11]);

    // Define local variables
    bool pass = false;

    // Launch test for zeroth order Legendre polynomial
    pass = launch_test(file_path_legendre0, 0);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 0 failed!" << std::endl;
        return 1;
    }

    // Launch test for first order Legendre polynomial
    pass = launch_test(file_path_legendre1, 1);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 1 failed!" << std::endl;
        return 1;
    }

    // Launch test for second order Legendre polynomial
    pass = launch_test(file_path_legendre2, 2);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 2 failed!" << std::endl;
        return 1;
    }

    // Launch test for third order Legendre polynomial
    pass = launch_test(file_path_legendre3, 3);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 3 failed!" << std::endl;
        return 1;
    }

    // Launch test for fourth order Legendre polynomial
    pass = launch_test(file_path_legendre4, 4);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 4 failed!" << std::endl;
        return 1;
    }

    // Launch test for fifth order Legendre polynomial
    pass = launch_test(file_path_legendre5, 5);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 4 failed!" << std::endl;
        return 1;
    }

    // Launch test for sixth order Legendre polynomial
    pass = launch_test(file_path_legendre6, 6);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 6 failed!" << std::endl;
        return 1;
    }

    // Launch test for seventh order Legendre polynomial
    pass = launch_test(file_path_legendre7, 7);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 7 failed!" << std::endl;
        return 1;
    }

    // Launch test for eight order Legendre polynomial
    pass = launch_test(file_path_legendre8, 8);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 8 failed!" << std::endl;
        return 1;
    }

    // Launch test for ninth order Legendre polynomial
    pass = launch_test(file_path_legendre9, 9);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 9 failed!" << std::endl;
        return 1;
    }

    // Launch test for tenth order Legendre polynomial
    pass = launch_test(file_path_legendre10, 10);
    if (!pass)
    {
        std::cerr << "test_legendre_poly/legendre_poly Order 10 failed!" << std::endl;
        return 1;
    }

    return 0;
}