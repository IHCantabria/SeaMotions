
// Include geneal usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math_tools.hpp" 


int sub_test_sv_add(cusfloat* u, cusfloat* v)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {5.0, 7.0, 9.0};

    // Perform vector addition
    sv_add(3, u, v, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_svs_add(cusfloat* u, cusfloat s)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {4.5, 5.5, 6.5};

    // Perform vector addition
    svs_add(3, u, s, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_sv_div(cusfloat* u, cusfloat* v)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {0.25, 0.4, 0.5};

    // Perform vector addition
    sv_div(3, u, v, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_sv_inv(cusfloat* u, cusfloat s)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {3.5, 1.75, 3.5/3};

    // Perform vector addition
    sv_inv(3, s, u, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_sv_mult(cusfloat* u, cusfloat* v)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {4.0, 10.0, 18.0};

    // Perform vector addition
    sv_mult(3, u, v, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_svs_mult(cusfloat* u, cusfloat s)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {3.5, 7.0, 10.5};

    // Perform vector addition
    svs_mult(3, u, s, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_sv_pow2(cusfloat* u)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {1.0, 4.0, 9.0};

    // Perform vector addition
    sv_pow2(3, u, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_sv_sqrt(cusfloat* u)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {1.0, std::sqrt(2.0), std::sqrt(3.0)};

    // Perform vector addition
    sv_sqrt(3, u, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_sv_sub(cusfloat* u, cusfloat* v)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {-3.0, -3.0, -3.0};

    // Perform vector addition
    sv_sub(3, u, v, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_svs_sub(cusfloat* u, cusfloat s)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {-2.5, -1.5, -0.5};

    // Perform vector addition
    svs_sub(3, u, s, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int sub_test_svs_sub2(cusfloat* u, cusfloat s)
{
    // Define result vector
    cusfloat w[3] = {0.0, 0.0, 0.0};
    cusfloat w_sol[3] = {-0.5, -1.5, -2.5};

    // Perform vector addition
    svs_sub(3, s, u, w);

    // Compare computed and solution vector
    int pass = assert_vector_equality(3, w, w_sol, EPS_PRECISION);

    return pass;
}


int main(void)
{
    // Declare pass variable to check tests
    int pass;

    // Define vectors and scalar to perform the tests
    cusfloat s = 3.5;
    cusfloat u[3] = {1.0, 2.0, 3.0};
    cusfloat v[3] = {4.0, 5.0, 6.0};

    // Test vector-vector element-wise addition
    pass = sub_test_sv_add(u, v);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_sv_add failed!" << std::endl;
        return 1;
    }

    // Test scalar-vector element-wise addition
    pass = sub_test_svs_add(u, s);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_svs_add failed!" << std::endl;
        return 1;
    }

    // Test vector-vector element-wise division
    pass = sub_test_sv_div(u, v);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_sv_div failed!" << std::endl;
        return 1;
    }

    // Test scalar-vector element-wise division
    pass = sub_test_sv_inv(u, s);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_svs_div failed!" << std::endl;
        return 1;
    }

    // Test vector-vector element-wise multiplication
    pass = sub_test_sv_mult(u, v);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_sv_mult failed!" << std::endl;
        return 1;
    }

    // Test scalar-vector element-wise multiplication
    pass = sub_test_svs_mult(u, s);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_svs_mult failed!" << std::endl;
        return 1;
    }

    // Test vector element-wise power of 2
    pass = sub_test_sv_pow2(u);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_svs_pow2 failed!" << std::endl;
        return 1;
    }

    // Test vector element-wise sqrt
    pass = sub_test_sv_sqrt(u);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_sv_sqrt failed!" << std::endl;
        return 1;
    }

    // Test vector-scalar element-wise substraction
    pass = sub_test_sv_sub(u, v);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_sv_sub failed!" << std::endl;
        return 1;
    }

    // Test scalar-vector element-wise substraction
    pass = sub_test_svs_sub2(v, s);
    if (pass == 0)
    {
        std::cerr << "test_vector_arithmetics/sub_test_svs_sub2 failed!" << std::endl;
        return 1;
    }

    return 0;
}