
// Include general usage libraries
#include <iostream>
#include <sstream>
#include <string>

// Include local modules
#include "../../src/math/math_tools.hpp"


template<typename T>
bool launch_test(void)
{
    // Compose genera error message
    std::stringstream err_message;
    err_message << "Sign test failed for type: " << typeid(T).name();

    // Check for positive numbers
    T np = 58420;
    if (sign(np) != 1)
    {
        err_message << " - Positive number: " << np << std::endl;
        std::cerr << err_message.str() << std::endl;
        return false;
    }

    // Check for zeros
    T nz = 0;
    if (sign(nz) != 0)
    {
        err_message << " - Zero number: " << np << std::endl;
        std::cerr << err_message.str() << std::endl;
        return false;
    }

    // Check for negative numbers
    T nn = -98541;
    if (sign(nn) != -1)
    {
        err_message << " - Negative number: " << np << std::endl;
        std::cerr << err_message.str() << std::endl;
        return false;
    }

    return true;
}


int main(void)
{
    // Declare pass flag
    bool pass = false;

    // Check for integer type
    pass = launch_test<int>();
    if (!pass)
    {
        return 1;
    }

    // Check for float type
    pass = launch_test<float>();
    if (!pass)
    {
        return 1;
    }

    // Check for double type
    pass = launch_test<double>();
    if (!pass)
    {
        return 1;
    }
    

    return 0;
}