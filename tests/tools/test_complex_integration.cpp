
// Include general usage libraries
#include <iostream>
#include <string>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"


// Include namespace
using namespace std::literals::complex_literals;


void assert_test(cuscomplex zn, cuscomplex za, std::string test_name, bool &pass)
{
    if (!assert_complex_equality(zn, za, EPS_PRECISION))
    {
        std::cout << "test_complex_integration/" << test_name << " failed!" << std::endl;
        std::cout << " - Int(num): " << zn << std::endl;
        std::cout << " - Int(ana): " << za << std::endl;
        std::cout << " - Diff(zn-za): " << zn-za << std::endl;
        pass = false;
    }
}


cuscomplex fcn_0(cuscomplex z)
{
    return z*z;
}


cuscomplex fcn_1(cuscomplex z)
{
    return 1.0/z;
}


cuscomplex fcn_2(cuscomplex z)
{
    return 1.0/((z*z)-1.0);
}


cuscomplex fcn_3(cuscomplex z)
{
    return 1.0/((z*z)+1.0);
}


bool launch_test0(void)
{
    // Define pass flag
    bool pass = true;

    // Define path of integration
    const int N = 2;
    cuscomplex way_points[N];
    way_points[0] = 0.0 + 0.0i;
    way_points[1] = 1.0 + 1.0i;

    // Integrate function along path
    cuscomplex zn = complex_integration_path(fcn_0, N, way_points, 1e-6, false);

    // Calculate theoretical solution
    cuscomplex ri = 1.0 + 1.0i;
    cuscomplex za = pow3s(ri)/3.0;

    // Compare numerical and analytical solutions
    assert_test(zn, za, "sub_test_0", pass);

    return pass;
}


bool launch_test1(void)
{
    // Declare test variable
    bool pass = false;

    // Define waypoints
    const int N = 17;
    cuscomplex waypoints[N];
    waypoints[0] = -1e6 + 0.0i;
    waypoints[1] = -1e5 + 0.0i;
    waypoints[2] = -1e4 + 0.0i;
    waypoints[3] = -1e3 + 0.0i;
    waypoints[4] = -1e2 + 0.0i;
    waypoints[5] = -1e1 + 0.0i;
    waypoints[6] = -2.0 + 0.0i;
    waypoints[7] = -1.0 + 0.0i;
    waypoints[8] =  0.0 + 1.0i;
    waypoints[9] =  1.0 + 0.0i;
    waypoints[10] = 2.0 + 0.0i;
    waypoints[11] = 1e1 + 0.0i;
    waypoints[12] = 1e2 + 0.0i;
    waypoints[13] = 1e3 + 0.0i;
    waypoints[14] = 1e4 + 0.0i;
    waypoints[15] = 1e5 + 0.0i;
    waypoints[16] = 1e6 + 0.0i;

    // Calculate Numerical Integral value
    cuscomplex zn = complex_integration_path(fcn_1, N, waypoints, 1e-12, false);

    // Set Analytical value
    cuscomplex za = 0.0 - PI*1i;

    // Compare Numerical and Analytical results
    assert_test(zn, za, "sub_test_1", pass);

    return pass;
}


bool launch_test2(void)
{
    // Declare test variable
    bool pass = false;

    // Define waypoints
    const int N = 17;
    cuscomplex waypoints[N];
    waypoints[0] = -1e6 + 0.0i;
    waypoints[1] = -1e5 + 0.0i;
    waypoints[2] = -1e4 + 0.0i;
    waypoints[3] = -1e3 + 0.0i;
    waypoints[4] = -1e2 + 0.0i;
    waypoints[5] = -1e1 + 0.0i;
    waypoints[6] = -1.5 + 0.0i;
    waypoints[7] = -1.0 + 1.0i;
    waypoints[8] = -0.5 + 0.0i;
    waypoints[8] =  0.5 + 0.0i;
    waypoints[9] =  1.0 + 1.0i;
    waypoints[10] =  1.5 + 0.0i;
    waypoints[11] =  1e1 + 0.0i;
    waypoints[12] =  1e2 + 0.0i;
    waypoints[13] =  1e3 + 0.0i;
    waypoints[14] =  1e4 + 0.0i;
    waypoints[15] =  1e5 + 0.0i;
    waypoints[16] =  1e6 + 0.0i;

    // Calulate Integral value
    cuscomplex z0 = complex_integration_path(fcn_2, N, waypoints, 1e-12, false);
    std::cout << "z0: " << z0 << std::endl;
    // cuscomplex a = 1.0 + 0.0i;
    // cuscomplex b = 2.0 + 0.0i;
    // cuscomplex z0 = complex_integration(fcn_1, a, b, 1e-6);
    // std::cout << "z0: " << z0 << std::endl;


    return pass;
}


bool launch_test3(void)
{
    // Declare test variable
    bool pass = false;

    // Define waypoints
    const int N = 19;
    cuscomplex waypoints[N];
    waypoints[0] = -1e7 + 0.0i;
    waypoints[1] = -1e5 + 0.0i;
    waypoints[2] = -1e4 + 0.0i;
    waypoints[3] = -1e3 + 0.0i;
    waypoints[4] = -1e2 + 0.0i;
    waypoints[5] = -1e1 + 0.0i;
    waypoints[6] = -2.0 + 0.0i;
    waypoints[7] = -1.0 + 0.0i;
    waypoints[8] =  0.0 + 0.0i;
    waypoints[9] =  1.0 + 0.0i;
    waypoints[10] =  2.0 + 0.0i;
    waypoints[11] =  1e1 + 0.0i;
    waypoints[12] =  1e2 + 0.0i;
    waypoints[13] =  1e3 + 0.0i;
    waypoints[14] =  1e4 + 0.0i;
    waypoints[15] =  1e5 + 0.0i;
    waypoints[16] =  1e7 + 0.0i;
    waypoints[17] =  0.0 + 5.0i;
    waypoints[18] = -1e7 + 0.0i;


    // Calulate Integral value
    cuscomplex z0 = complex_integration_path(fcn_2, N, waypoints, 1e-6, false);
    std::cout << "z0: " << z0 << std::endl;
    // cuscomplex a = 1.0 + 0.0i;
    // cuscomplex b = 2.0 + 0.0i;
    // cuscomplex z0 = complex_integration(fcn_1, a, b, 1e-6);
    // std::cout << "z0: " << z0 << std::endl;


    return pass;
}


int main(void)
{
    // Declare logic system variables
    bool pass = false;

    // Launch test 0 -> int(z^2, 0.0+0.0i, 1.0+1.0i)
    // pass = launch_test0();
    // if (!pass)
    // {
    //     return 1;
    // }

    // // Launch test 1 -> int(1/z, -inf+0.0i, inf+0.0i)
    // pass = launch_test1();
    // if (!pass)
    // {
    //     return 1;
    // }

    // Launch test 2 -> int(1/(z^2-1), -inf+0.0i, inf+0.0i)
    pass = launch_test2();
    if (!pass)
    {
        return 1;
    }

    // // Launch test 3 -> int(1/(z^2+1), -inf+0.0i, inf+0.0i)
    // pass = launch_test3();
    // if (!pass)
    // {
    //     return 1;
    // }

    return 0;
}