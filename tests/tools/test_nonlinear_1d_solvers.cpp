
// Include general usage libraries
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/math/math_tools.hpp"

// Define tests precision
cusfloat ABS_TOL = 1e-6;
cusfloat REL_TOL = 0.0;


cusfloat fcn_1(cusfloat x)
{
    return pow2s(x)-1.0;
}


cusfloat fcn_1_der(cusfloat x)
{
    return 2*x;
}


cusfloat fcn_2(cusfloat x)
{
    return std::sin(x)-x/4;
}


cusfloat fcn_2_der(cusfloat x)
{
    return std::cos(x)-1/4;
}


cusfloat fcn_3(cusfloat x)
{
    return std::atanh(x)-1.0;
}


cusfloat fcn_3_der(cusfloat x)
{
    return 1/(1-pow2s(x));
}


bool launch_bisection(std::function<cusfloat(cusfloat)> f_def, cusfloat a, cusfloat b, cusfloat xt)
{
    bool pass = true;
    int info = -1;
    cusfloat xs = 0.0;
    bisection(f_def, a, b, ABS_TOL, REL_TOL, 100, false, xs, info);
    if ((abs(xt-xs)>ABS_TOL) || (info != 0))
    {
        pass = false;
    }

    return pass;
}


bool launch_newton(std::function<cusfloat(cusfloat)> f_def, std::function<cusfloat(cusfloat)> f_der_def,
    cusfloat x0, cusfloat xt)
{
    bool pass = true;
    int info = -1;
    cusfloat xs = 0.0;
    newton_raphson(f_def, f_der_def, x0, ABS_TOL, REL_TOL, 100, false, xs, info);
    if ((abs(xt-xs)>ABS_TOL) || (info != 0))
    {
        pass = false;
    }

    return pass;
}


bool bisection_tests(void)
{
    // Define pass flag
    bool pass = true;

    // Define common error string name
    std::string err_common("test_nonlinear_solvers/bisection_tests/");
    std::string err_common_fin(" failed!");

    // Find solution using bisection method to fcn_1
    cusfloat xt = 1.0;
    pass = launch_bisection(fcn_1, 0.01, 1.84, xt);
    if (!pass)
    {
        std::cerr << err_common << "fcn_1" << err_common_fin << std::endl;
        return pass;
    }

    // Find solution using bisection method to fcn_2
    xt = 2.474576787369829;
    pass = launch_bisection(fcn_2, 2.0, 3.0, xt);
    if (!pass)
    {
        std::cerr << err_common << "fcn_2" << err_common_fin << std::endl;
        return pass;
    }

    // Find solution using bisection method to fcn_2
    xt = 0.761594155955765;
    pass = launch_bisection(fcn_3, 0.5, 0.92, xt);
    if (!pass)
    {
        std::cerr << err_common << "fcn_3" << err_common_fin << std::endl;
        return pass;
    }

    return pass;
}

bool newton_tests(void)
{
    // Define pass flag
    bool pass = true;

    // Define common error string name
    std::string err_common("test_nonlinear_solvers/newton_tests/");
    std::string err_common_fin(" failed!");

    // Find solution using bisection method to fcn_1
    cusfloat xt = 1.0;
    pass = launch_newton(fcn_1, fcn_1_der, 0.01, xt);
    if (!pass)
    {
        std::cerr << err_common << "fcn_1" << err_common_fin << std::endl;
        return pass;
    }

    // Find solution using bisection method to fcn_2
    xt = 2.474576787369829;
    pass = launch_newton(fcn_2, fcn_2_der, 2.0, xt);
    if (!pass)
    {
        std::cerr << err_common << "fcn_2" << err_common_fin << std::endl;
        return pass;
    }

    // Find solution using bisection method to fcn_2
    xt = 0.761594155955765;
    pass = launch_newton(fcn_3, fcn_3_der, 0.5, xt);
    if (!pass)
    {
        std::cerr << err_common << "fcn_3" << err_common_fin << std::endl;
        return pass;
    }

    return pass;
}


int main(void)
{
    // Declare system logic flags
    int pass = false;
    int return_chn = 0;

    // Launch tests for bisection method
    pass = bisection_tests();
    if (!pass)
    {
        return_chn = 1;
    }

    // Launch tests for bisection method
    pass = newton_tests();
    if (!pass)
    {
        return_chn = 1;
    }


    return return_chn;
}