
// Include external libraries
#include <cmath>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/tools.hpp"


void print_vector(cusfloat (&v)[3])
{
    std::cout << "v[" << 0 << "]: " << v[0];
    for (int i=1; i<3; i++)
    {
        std::cout << " - v[" << i << "]: " << v[i];
    }
    std::cout << std::endl;

}


int sub_test_1(void)
{
    // Define test vectors
    cusfloat v0[3] = {1.0, 0.0, 0.0};
    cusfloat v1[3] = {0.0, 1.0, 0.0};
    cusfloat v2[3] = {0.0, 0.0, 0.0};

    // Define solution  vector
    cusfloat v3[3] = {0.0, 0.0, 1.0};

    // Excecute cross prouct
    cross(v0, v1, v2);

    // Check if the calculated vector is equal to the expected
    int result = assert_vector_equality<cusfloat>(v2, v3, EPS_PRECISION);

    return result;
}


int sub_test_2(void)
{
    // Define test vectors
    cusfloat v0[3] = {0.0, 1.0, 0.0};
    cusfloat v1[3] = {0.0, 0.0, 1.0};
    cusfloat v2[3] = {0.0, 0.0, 0.0};

    // Define solution  vector
    cusfloat v3[3] = {1.0, 0.0, 0.0};

    // Excecute cross prouct
    cross(v0, v1, v2);

    // Check if the calculated vector is equal to the expected
    int result = assert_vector_equality<cusfloat>(v2, v3, EPS_PRECISION);

    return result;
}


int sub_test_3(void)
{
    // Define test vectors
    cusfloat v0[3] = {0.0, 0.0, 1.0};
    cusfloat v1[3] = {1.0, 0.0, 0.0};
    cusfloat v2[3] = {0.0, 0.0, 0.0};

    // Define solution  vector
    cusfloat v3[3] = {0.0, 1.0, 0.0};

    // Excecute cross prouct
    cross(v0, v1, v2);

    // Check if the calculated vector is equal to the expected
    int result = assert_vector_equality<cusfloat>(v2, v3, EPS_PRECISION);

    return result;
}


int sub_test_4(void)
{
    // Define test vectors
    int v0[3] = {0, 0, 1};
    int v1[3] = {1, 0, 0};
    int v2[3] = {0, 0, 0};

    // Define solution  vector
    int v3[3] = {0, 1, 0};

    // Excecute cross prouct
    cross(v0, v1, v2);

    // Check if the calculated vector is equal to the expected
    int result = assert_vector_equality<int>(v2, v3, 0);

    return result;
}


int sub_test_5(void)
{
    // Define test vectors
    cusfloat v0[3] = {0.250, 0.368, -0.980};
    cusfloat v1[3] = {-1.1, -0.0032, 0.7};
    cusfloat v2[3] = {0, 0, 0};

    // Define solution  vector
    cusfloat v3[3] = {0.254464, 0.903, 0.404};

    // Excecute cross prouct
    cross(v0, v1, v2);

    // Check if the calculated vector is equal to the expected
    int result = assert_vector_equality<cusfloat>(v2, v3, EPS_PRECISION);

    return result;
}


int main(void)
{
    //Define and launch first sub-test
    int result_1 = sub_test_1();
    if (result_1 == 0)
    {
        std::cerr << "test_dipole_potentail/sub_test_1 failed!" << std::endl;
        return 1;
    }
    
    //Define and launch second sub-test
    int result_2 = sub_test_2();
    if (result_2 == 0)
    {
        std::cerr << "test_dipole_potentail/sub_test_2 failed!" << std::endl;
        return 1;
    }

    //Define and launch thrid sub-test
    int result_3 = sub_test_3();
    if (result_3 == 0)
    {
        std::cerr << "test_dipole_potentail/sub_test_3 failed!" << std::endl;
        return 1;
    }

    //Define and launch fourth sub-test
    int result_4 = sub_test_4();
    if (result_4 == 0)
    {
        std::cerr << "test_dipole_potentail/sub_test_4 failed!" << std::endl;
        return 1;
    }

    //Define and launch fifth sub-test
    int result_5 = sub_test_5();
    if (result_5 == 0)
    {
        std::cerr << "test_dipole_potentail/sub_test_5 failed!" << std::endl;
        return 1;
    }

    return 0;
}