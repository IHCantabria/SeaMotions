
// Include general usage libraries
#include <functional>
#include <iostream>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating.hpp"
#include "../../src/math_tools.hpp"
#include "../../src/special_math.hpp"


// Declare local module functions
void generate_grid(int n, cusfloat y_min, cusfloat y_max, cusfloat* x, cusfloat* y);
void get_x_range(cusfloat y, cusfloat &x0, cusfloat &x1);


cusfloat eval_chebyshev_xy(int n_order, cusfloat x, cusfloat y, cusfloat* cheby_coeff)
{
    cusfloat sol = 0.0;
    for (int i=0; i<n_order; i++)
    {
        for (int j=0; j<n_order; j++)
        {
            sol += cheby_coeff[i*n_order+j]*chebyshev_poly_raw(i, x)*chebyshev_poly_raw(j, y);
            // std::cout << "i: " << i << " - j: " << j << " - sol: " << sol << " - chev: " << cheby_coeff[i*n_order+j] << std::endl;
        }
    }

    return sol;
}


void fit_chebyshev_xy(int n_order, int num_points, cusfloat* x, cusfloat * y,
    cusfloat* f, cusfloat* cheby_coeff)
{
    // Generate kernel matrix
    int num_cols = n_order*n_order;
    cusfloat* A = generate_empty_vector<cusfloat>(num_points*num_cols);
    cusfloat* At = generate_empty_vector<cusfloat>(num_points*num_cols);
    for (int i=0; i<num_points; i++)
    {
        for (int j=0; j<n_order; j++)
        {
            for (int k=0; k<n_order; k++)
            {
                A[i*num_cols+j*n_order+k] = chebyshev_poly_raw(j, x[i])*chebyshev_poly_raw(k, y[i]);
                At[(j*n_order+k)*num_points+i] = A[i*num_cols+j*n_order+k];
            }
        }
    }

    cusfloat* Ak = generate_empty_vector<cusfloat>(num_cols*num_cols);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                num_cols, num_cols, num_points, 1.0, At, num_points, A, num_cols, 1.0, Ak, num_cols);

    // Generate projected target vector
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                num_cols, num_points, 1.0, At, num_points, f, 1, 0.0, cheby_coeff, 1);

    // Solve linear system of equations to get the Chebyshev series coefficients
    int* ipiv = generate_empty_vector<int>(num_cols);
    int info = 0;
    int num_rhs = 1;
    dgesv(&num_cols, &num_rhs, Ak, &num_cols, ipiv, cheby_coeff, &num_cols, &info);

    // Delete heap memory
    mkl_free(A);
    mkl_free(Ak);
    mkl_free(At);
    mkl_free(ipiv);
}


void generate_grid(int n, cusfloat y_min, cusfloat y_max, cusfloat* x, cusfloat* y)
{
    // Generate y coordinate grid
    cusfloat dy = (y_max-y_min)/(n-1);
    cusfloat yi = 0.0;
    for (int i=0; i<n; i++)
    {
        yi = y_min + i*dy;
        for (int j=0; j<n; j++)
        {
            y[i*n+j] = yi;
        }
    }

    // Generate x coordinate grid
    cusfloat dx = 0.0;
    cusfloat x0 = 0.0, x1 = 0.0;

    for (int i=0; i<n; i++)
    {
        get_x_range(y[i*n], x0, x1);
        dx = (x1-x0)/(n-1);
        for (int j=0; j<n; j++)
        {
            x[i*n+j] = x0 + j*dx;
        }
    }

}


void get_x_range(cusfloat y, cusfloat &x0, cusfloat &x1)
{
    x0 = y/2.0;
    x1 = 4.0*y;
}


void fit_sub_domain_a(void)
{
    // Define allocation matrixes
    constexpr int num_points = 50;
    constexpr int num_points_p2 = num_points*num_points;
    cusfloat x[num_points_p2], y[num_points_p2];

    // Generate grid
    generate_grid(num_points, 2, 4, x, y);

    // std::cout << "X Coordinate" << std::endl;
    // print_matrix(num_points, num_points, x, 3, 0, 0);

    // std::cout << "Y Coordinate" << std::endl;
    // print_matrix(num_points, num_points, y, 3, 0, 0);

    // Generate target vector
    cusfloat* f = generate_empty_vector<cusfloat>(num_points);
    for (int i=0; i<num_points; i++)
    {
        f[i] = wave_term_inf_depth_num(x[i], y[i]);
    }

    // Fit Chebyshev polynomials
    constexpr int n_order = 5;
    cusfloat cheby_coeff[n_order*n_order];
    fit_chebyshev_xy(n_order, num_points_p2, x, y, f, cheby_coeff);

    for (int i=0; i<n_order*n_order; i++)
    {
        std::cout << "C[" << i << "]:" << cheby_coeff[i] << std::endl;
    }

    // Evaluate Chebyshev polynomials
    cusfloat* sol_fit = generate_empty_vector<cusfloat>(num_points_p2);
    for (int i=0; i<=num_points_p2; i++)
    {
        sol_fit[i] = eval_chebyshev_xy(n_order, x[i], y[i], cheby_coeff);
    }

    // Compute error with respecto the reference values
    cusfloat* err = generate_empty_vector<cusfloat>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        err[i] = sol_fit[i]-f[i];
    }

    for (int i=0; i<num_points_p2; i++)
    {
        std::cout << "X: " << x[i] << " - Y: " << y[i];
        std::cout << " - fi: " << f[i] << " - sol_fit[i]: " << sol_fit[i];
        std::cout << " - err[" << i << "]: " << err[i] << std::endl;
    }

    // Delete heap memory
    mkl_free(f);
    mkl_free(sol_fit);

}


int main(void)
{
    fit_sub_domain_a();

    return 0;
}