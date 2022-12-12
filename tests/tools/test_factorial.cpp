
// Include general usage libraries
#include <iostream>

// Include local modules
#include "../../src/math_tools.hpp"


int sub_test_1(void)
{
    int pass = 1;
    if (std::abs(1-factorial(0)) > 0.0)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_2(void)
{
    int pass = 1;
    if (std::abs(1-factorial(1)) > 0.0)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_3(void)
{
    int pass = 1;
    if (std::abs(2-factorial(2)) > 0.0)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_4(void)
{
    int pass = 1;
    if (std::abs(5040-factorial(7)) > 0.0)
    {
        pass = 0;
    }

    return pass;
}


int sub_test_5(void)
{
    int pass = 1;
    if (std::abs(479001600-factorial(12)) > 0.0)
    {
        pass = 0;
    }

    return pass;
}


int main(void)
{
    // Declare local variables
    int pass = 0;

    // Define first test
    pass = sub_test_1();
    if (pass == 0)
    {
        std::cerr << "test_factorial/sub_test_1 failed!" << std::endl;
        return 1;
    }

    // Define second test
    pass = sub_test_2();
    if (pass == 0)
    {
        std::cerr << "test_factorial/sub_test_2 failed!" << std::endl;
        return 1;
    }

    // Define third test
    pass = sub_test_3();
    if (pass == 0)
    {
        std::cerr << "test_factorial/sub_test_3 failed!" << std::endl;
        return 1;
    }

    // Define fourth test
    pass = sub_test_4();
    if (pass == 0)
    {
        std::cerr << "test_factorial/sub_test_4 failed!" << std::endl;
        return 1;
    }

    // Define fifth test
    pass = sub_test_5();
    if (pass == 0)
    {
        std::cerr << "test_factorial/sub_test_5 failed!" << std::endl;
        return 1;
    }


    return 0;
}