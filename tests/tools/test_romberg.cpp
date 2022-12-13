
// Include external libraries
#include <iostream>

// Include external scientific libraries
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math_tools.hpp"


cusfloat NUM_EPS_PRECISION = 1e-8;


cusfloat bessel_fcn(cusfloat x)
{
    return -1*_j1(x);
}


cusfloat linear_fcn(cusfloat x)
{
    return x;
}


cusfloat poly_fcn(cusfloat x)
{
    return 3*pow3s(x)+4*pow2s(x)+x+5.0;
}


cusfloat poly_int_fcn(cusfloat x)
{
    return 3.0/4.0*std::pow(x, 4)+4.0/3.0*pow3s(x)+pow2s(x)/2.0+5.0*x;
}


cusfloat sqrt_fcn(cusfloat x)
{
    return 0.5/std::sqrt(x);
}


int sub_test_1(void)
{
    // Define analytical solution
    cusfloat ana_sol = 0.5;

    // Calculate numerical solution using Romberg method
    cusfloat num_sol = romberg_quadrature(linear_fcn, 0, 1, NUM_EPS_PRECISION);

    // Check the accurcy of the solution
    int pass = 1;
    if (std::abs(ana_sol-num_sol) > NUM_EPS_PRECISION)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_2(void)
{
    // Define analytical solution
    cusfloat a = -7.1932;
    cusfloat b = 5.418;
    cusfloat ana_sol = poly_int_fcn(b) - poly_int_fcn(a);

    // Calculate numerical solution using Romberg method
    cusfloat num_sol = romberg_quadrature(poly_fcn, a, b, NUM_EPS_PRECISION);

    // Check the accurcy of the solution
    int pass = 1;
    if (std::abs(ana_sol-num_sol) > NUM_EPS_PRECISION)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_3(void)
{
    // Define analytical solution
    cusfloat ana_sol = std::exp(1.0) - std::exp(0.0);

    // Calculate numerical solution using Romberg method
    cusfloat num_sol = romberg_quadrature([](cusfloat t) -> cusfloat {return std::exp(t);}, 0.0, 1.0, 1e-15);

    // Check the accurcy of the solution
    int pass = 1;
    if (std::abs(ana_sol-num_sol) > NUM_EPS_PRECISION)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_4(void)
{
    // Define analytical solution
    cusfloat a = 0.01425;
    cusfloat b = 6.5214;
    cusfloat ana_sol = std::sqrt(b) - std::sqrt(a);

    // Calculate numerical solution using Romberg method
    cusfloat num_sol = romberg_quadrature(sqrt_fcn, a, b, NUM_EPS_PRECISION);

    // Check the accurcy of the solution
    int pass = 1;
    if (std::abs(ana_sol-num_sol) > NUM_EPS_PRECISION)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_5(void)
{
    // Define analytical solution
    cusfloat a = 0.01425;
    cusfloat b = 6.5214;
    cusfloat ana_sol = _j0(b) - _j0(a);

    // Calculate numerical solution using Romberg method
    cusfloat num_sol = romberg_quadrature(bessel_fcn, a, b, NUM_EPS_PRECISION);

    // Check the accurcy of the solution
    int pass = 1;
    if (std::abs(ana_sol-num_sol) > NUM_EPS_PRECISION)
    {
        pass = 0;
    }

    return pass;
}


int main(void)
{
    int pass;

    // Check romberg accuracy with simple linear function
    pass = sub_test_1();
    if (pass == 0)
    {
        std::cerr << "test_romberg/sub_test_1 failed!" << std::endl;
        return 1;
    }

    // Check romberg accuracy with a polynomial function
    pass = sub_test_2();
    if (pass == 0)
    {
        std::cerr << "test_romberg/sub_test_2 failed!" << std::endl;
        return 1;
    }

    // Check romberg accuracy with an exponential function
    pass = sub_test_3();
    if (pass == 0)
    {
        std::cerr << "test_romberg/sub_test_3 failed!" << std::endl;
        return 1;
    }

    // Check romberg accuracy with inverse sqrt function
    pass = sub_test_4();
    if (pass == 0)
    {
        std::cerr << "test_romberg/sub_test_4 failed!" << std::endl;
        return 1;
    }

    // Check romberg accuracy with the bessel function of the first kind and first order
    pass = sub_test_5();
    if (pass == 0)
    {
        std::cerr << "test_romberg/sub_test_5 failed!" << std::endl;
        return 1;
    }

    return 0;
}