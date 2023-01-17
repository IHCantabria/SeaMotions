
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/math/math_tools.hpp"

// Define tests precision
cusfloat ABS_PREC = 1e-6;


cusfloat fcn_1(cusfloat x)
{
    return pow2s(x)-1.0;
}


cusfloat fcn_2(cusfloat x)
{
    return std::sin(x)-x/4;
}


cusfloat fcn_3(cusfloat x)
{
    return std::atanh(x)-1.0;
}


bool sub_test_1(void)
{
    // Define pass flag
    bool pass = true;

    // Find solution using bisection method to fcn_1
    cusfloat xt = 1.0;
    cusfloat xs = bisection(fcn_1, 0.01, 1.84, ABS_PREC, 100, false);
    if (abs(xt-xs)>ABS_PREC)
    {
        pass = false;
    }

    return pass;
}


bool sub_test_2(void)
{
    // Define pass flag
    bool pass = true;

    // Find solution using bisection method to fcn_1
    cusfloat xt = 2.474576787369829;
    cusfloat xs = bisection(fcn_2, 2.0, 3.0, ABS_PREC, 100, false);
    if (abs(xt-xs)>ABS_PREC)
    {
        pass = false;
    }

    return pass;
}


bool sub_test_3(void)
{
    // Define pass flag
    bool pass = true;

    // Find solution using bisection method to fcn_1
    cusfloat xt = 0.761594155955765;
    cusfloat xs = bisection(fcn_3, 0.5, 0.92, ABS_PREC, 100, false);
    if (abs(xt-xs)>ABS_PREC)
    {
        pass = false;
    }

    return pass;
}


int main(void)
{
    // Declare pass flat
    int pass = false;

    // Launch test for: f(x)=x^2-1=0
    pass = sub_test_1();
    if (!pass)
    {
        std::cerr << "test_bisection/sub_test_1 failed!" << std::endl;
    }

    // Launch test for: f(x)=sin(x)-x/4=0
    pass = sub_test_2();
    if (!pass)
    {
        std::cerr << "test_bisection/sub_test_2 failed!" << std::endl;
    }

    // Launch test for: f(x)=atan(x)-1=0
    pass = sub_test_3();
    if (!pass)
    {
        std::cerr << "test_bisection/sub_test_3 failed!" << std::endl;
    }


    return 0;
}