
// Include general usage libraries
#include <iostream>


int sub_test_1(void)
{
    

    return 0;
}


int main(void)
{
    int pass;

    // Compare velocity field calculated through Hess and Smith with Newman formulation
    pass = sub_test_1();
    if (pass == 0)
    {
        std::cerr << "test_source_velocity/sub_test_1 failed!" << std::endl;
        return 1;
    }

    return 0;
}