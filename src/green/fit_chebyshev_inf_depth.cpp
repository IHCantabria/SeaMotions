
// Include general usage libraries
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating.hpp"
#include "../../src/math_tools.hpp"
#include "../../src/special_math.hpp"

// Redefine-include name spaces
namespace fs = std::filesystem;


// Declare local module functions
void generate_grid(int n, cusfloat y_min, cusfloat y_max, cusfloat* x, cusfloat* y);
void get_x_range(cusfloat y, cusfloat &x0, cusfloat &x1, int side);


struct FitRegion
{
    cusfloat dy = 0.0;
    int num_points = 0;
    int side = 0;
    cusfloat y_min = 0.0;
    cusfloat y_max = 0.0;

    void fit(cusfloat x, cusfloat y, cusfloat &xl, cusfloat &yl)
    {
        // Calculate x local range
        cusfloat x0 = 0.0, x1 = 0.0;
        get_x_range(y, x0, x1, this->side);
        cusfloat dx = x1-x0;
        xl = 2.0*(x-x0)/dx-1.0;

        // Calculate y local range
        yl = 2.0*(y-this->y_min)/this->dy-1.0;

    }

};


cusfloat eval_chebyshev_xy(int n_order, FitRegion fr, cusfloat x, cusfloat y, cusfloat* cheby_coeff)
{
    cusfloat sol = 0.0;
    cusfloat xl = 0.0, yl = 0.0;
    for (int i=0; i<n_order; i++)
    {
        for (int j=0; j<n_order; j++)
        {
            fr.fit(x, y, xl, yl);
            sol += cheby_coeff[i*n_order+j]*chebyshev_poly_raw(i, xl)*chebyshev_poly_raw(j, yl);
            // std::cout << "i: " << i << " - j: " << j << " - sol: " << sol << " - chev: " << cheby_coeff[i*n_order+j] << std::endl;
        }
    }

    return sol;
}


void fit_chebyshev_xy(int n_order, int num_points, FitRegion fr, cusfloat* x, cusfloat * y,
    cusfloat* f, cusfloat* cheby_coeff)
{
    // Generate kernel matrix
    int num_cols = n_order*n_order;
    cusfloat* A = generate_empty_vector<cusfloat>(num_points*num_cols);
    cusfloat* At = generate_empty_vector<cusfloat>(num_points*num_cols);
    cusfloat xl = 0.0, yl = 0.0;
    for (int i=0; i<num_points; i++)
    {
        fr.fit(x[i], y[i], xl, yl);
        // std::cout << "X: " << x[i] << " - Y[i]: " << y[i] << " - xl: " << xl << " - yl: " << yl << std::endl;
        for (int j=0; j<n_order; j++)
        {
            for (int k=0; k<n_order; k++)
            {
                A[i*num_cols+j*n_order+k] = chebyshev_poly_raw(j, xl)*chebyshev_poly_raw(k, yl);
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


void generate_grid(FitRegion fr, cusfloat* x, cusfloat* y)
{
    // Generate y coordinate grid
    cusfloat dy = (fr.y_max-fr.y_min)/(fr.num_points-1);
    cusfloat yi = 0.0;
    for (int i=0; i<fr.num_points; i++)
    {
        yi = fr.y_min + i*dy;
        for (int j=0; j<fr.num_points; j++)
        {
            y[i*fr.num_points+j] = yi;
        }
    }

    // Generate x coordinate grid
    cusfloat dx = 0.0;
    cusfloat x0 = 0.0, x1 = 0.0;

    for (int i=0; i<fr.num_points; i++)
    {
        get_x_range(y[i*fr.num_points], x0, x1, fr.side);
        dx = (x1-x0)/(fr.num_points-1);
        for (int j=0; j<fr.num_points; j++)
        {
            x[i*fr.num_points+j] = x0 + j*dx;
        }
    }

}


void get_x_range(cusfloat y, cusfloat &x0, cusfloat &x1, int side)
{
    if (y<=4.0)
    {
        if (side == 0)
        {
            x0 = y/2.0;
            x1 = 3.0;
        }
        else
        {
            x0 = 3.0;
            x1 = 4.0*y;
        }
    }
    else if (y<=8.0)
    {
        if (side == 0)
        {
            x0 = y/2.0;
            x1 = 8.0;
        }
        else
        {
            x0 = 8.0;
            x1 = 4.0*y;
        }
    }
    else
    {
        x0 = y/2.0;
        x1 = 4.0*y;
    }
}


void fit_sub_domain(FitRegion fr, int cheby_order, int &count_coeffs,
    cusfloat* cheby_coeff_filter, int* cheby_order_0, 
    int* cheby_order_1, std::string sub_domain_name)
{
    // Define allocation matrixes
    int cheby_order_2 = cheby_order*cheby_order;
    int num_points_p2 = fr.num_points*fr.num_points;
    std::cout << "num_points_p2: " << num_points_p2 << std::endl;
    cusfloat* x = generate_empty_vector<cusfloat>(num_points_p2);
    cusfloat* y = generate_empty_vector<cusfloat>(num_points_p2);

    // Generate grid
    std::cout << "Generate grid... " << std::endl;
    generate_grid(fr, x, y);

    // Generate target vector
    std::cout << "Generate target vector... " << std::endl;
    cusfloat* f = generate_empty_vector<cusfloat>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        f[i] = (
            expint_inf_depth_num(x[i], y[i])
            -1/std::sqrt(pow2s(x[i])+pow2s(y[i]))
            +std::exp(-y[i])/x[i]
            )*pow3s(std::sqrt(pow2s(x[i])+pow2s(y[i])))/y[i];
    }

    cusfloat* fxy = generate_empty_vector<cusfloat>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        fxy[i] = wave_term_inf_depth_num(x[i], y[i]) + 1/std::sqrt(pow2s(x[i])+pow2s(y[i]));
    }

    // Fit Chebyshev polynomials
    cusfloat* cheby_coeff = generate_empty_vector<cusfloat>(cheby_order_2);
    std::cout << "Fit " << sub_domain_name << std::endl;
    fit_chebyshev_xy(cheby_order, num_points_p2, fr, x, y, f, cheby_coeff);

    // Evaluate Chebyshev polynomials
    cusfloat* sol_fit = generate_empty_vector<cusfloat>(num_points_p2);
    cusfloat f1 = 0.0;
    cusfloat R = 0.0;
    for (int i=0; i<=num_points_p2; i++)
    {
        R = std::sqrt(pow2s(x[i])+pow2s(y[i]));
        f1 = 1/R - std::exp(-y[i])/x[i] + y[i]/pow3s(R)*eval_chebyshev_xy(cheby_order, fr, x[i], y[i], cheby_coeff);
        f1 = 1/R - PI*std::exp(-y[i])*(struve0(x[i])+bessely0(x[i])) - 2.0*f1;
        sol_fit[i] = f1;
    }

    // Compute error with respecto the reference values
    cusfloat* err = generate_empty_vector<cusfloat>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        err[i] = sol_fit[i]-fxy[i];
    }

    // Compute error statistics
    int count_thr = 0;
    cusfloat max_err = 0.0;
    cusfloat mean_err = 0.0;
    cusfloat min_err = 0.0;
    cusfloat threshold = 1e-6;
    for (int i=0; i<num_points_p2; i++)
    {
        // Check for threshold value
        if (std::abs(err[i])>threshold)
        {
            count_thr++;
        }

        // Check for the maximum value
        if (std::abs(err[i])>std::abs(max_err))
        {
            max_err = err[i];
        }

        if (std::abs(err[i])<std::abs(min_err))
        {
            min_err = err[i];
        }

        mean_err += std::abs(err[i]);
    }
    mean_err /= num_points_p2;

    std::cout << sub_domain_name << " - Error statistics:" << std::endl;
    std::cout << " - Max: " << max_err << std::endl;
    std::cout << " - Min: " << min_err << std::endl;
    std::cout << " - Mean: " << mean_err << std::endl;
    std::cout << " - Over threshold [" << threshold << "]: " << count_thr*1.0/num_points_p2 << std::endl;
    std::cout << std::endl;

    // Filter Chebyshev coefficients by threshold
    cusfloat cheby_thres = 1e-7;
    for (int i=0; i<cheby_order_2; i++)
    {
        if (std::abs(cheby_coeff[i])>=cheby_thres)
        {
            count_coeffs++;
        }
    }

    cheby_coeff_filter = generate_empty_vector<cusfloat>(count_coeffs);
    cheby_order_0 = generate_empty_vector<int>(count_coeffs);
    cheby_order_1 = generate_empty_vector<int>(count_coeffs);
    int count_index = 0;
    for (int i=0; i<cheby_order; i++)
    {
        for (int j=0; j<cheby_order; j++)
        {
            if (std::abs(cheby_coeff[i*cheby_order+j])>cheby_thres)
            {
                cheby_coeff_filter[count_index] = cheby_coeff[i*cheby_order+j];
                cheby_order_0[count_index] = i;
                cheby_order_1[count_index] = j;
                count_index++;
            }
        }
    }
 
    // Delete heap memory
    mkl_free(cheby_coeff);
    mkl_free(f);
    mkl_free(fxy);
    mkl_free(sol_fit);
    mkl_free(err);

}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    fs::path folder_path(argv[1]);
    fs::path file_path = folder_path / "my_file.dat";

    std::cout << "Folder path: " << folder_path << std::endl;
    std::cout << "Folder path: " << file_path << std::endl;

    // Allocate heap memory to storage fit results
    cusfloat* cheby_coeff_filter_a0 = nullptr;
    cusfloat* cheby_coeff_filter_a1 = nullptr;
    cusfloat* cheby_coeff_filter_b0 = nullptr;
    cusfloat* cheby_coeff_filter_b1 = nullptr;
    cusfloat* cheby_coeff_filter_c = nullptr;
    int* cheby_order_0_a0 = nullptr;
    int* cheby_order_0_a1 = nullptr;
    int* cheby_order_0_b0 = nullptr;
    int* cheby_order_0_b1 = nullptr;
    int* cheby_order_0_c = nullptr;
    int* cheby_order_1_a0 = nullptr;
    int* cheby_order_1_a1 = nullptr;
    int* cheby_order_1_b0 = nullptr;
    int* cheby_order_1_b1 = nullptr;
    int* cheby_order_1_c = nullptr;
    int count_cheby_a0 = 0;
    int count_cheby_a1 = 0;
    int count_cheby_b0 = 0;
    int count_cheby_b1 = 0;
    int count_cheby_c = 0;

    // Define fit properties
    int num_points = 100;
    int cheby_order = 10;

    // Generate Chebyshev polynomials for region a0 (X<3.0 & 2.0<=Y<=4.0)
    std::string sub_domain_name_a0 = "Sub-Domain A0";
    FitRegion fr_a0;
    fr_a0.num_points = num_points;
    fr_a0.side = 0;
    fr_a0.y_min = 2.0;
    fr_a0.y_max = 4.0;
    fr_a0.dy = (fr_a0.y_max-fr_a0.y_min);
    fit_sub_domain(fr_a0, cheby_order, count_cheby_a0, cheby_coeff_filter_a0,
        cheby_order_0_a0, cheby_order_1_a0, sub_domain_name_a0);

    // Generate Chebyshev polynomials for region a1 (X>3.0 & 2.0<=Y<=4.0)
    std::string sub_domain_name_a1 = "Sub-Domain A1";
    FitRegion fr_a1;
    fr_a1.num_points = num_points;
    fr_a1.side = 1;
    fr_a1.y_min = 2.0;
    fr_a1.y_max = 4.0;
    fr_a1.dy = (fr_a1.y_max-fr_a1.y_min);
    fit_sub_domain(fr_a1, cheby_order, count_cheby_a1, cheby_coeff_filter_a1,
        cheby_order_0_a1, cheby_order_1_a1, sub_domain_name_a1);

    // Generate Chebyshev polynomials for region b0 (X<8.0 & 4.0<Y<=8.0)
    std::string sub_domain_name_b0 = "Sub-Domain B0";
    FitRegion fr_b0;
    fr_b0.num_points = num_points;
    fr_b0.side = 0;
    fr_b0.y_min = 4.0001;
    fr_b0.y_max = 8.0;
    fr_b0.dy = (fr_b0.y_max-fr_b0.y_min);
    fit_sub_domain(fr_b0, cheby_order, count_cheby_b0, cheby_coeff_filter_b0,
        cheby_order_0_b0, cheby_order_1_b0, sub_domain_name_b0);

    // Generate Chebyshev polynomials for region b1 (X>8.0 & 4.0<Y<=8.0)
    std::string sub_domain_name_b1 = "Sub-Domain B1";
    FitRegion fr_b1;
    fr_b1.num_points = num_points;
    fr_b1.side = 1;
    fr_b1.y_min = 4.0001;
    fr_b1.y_max = 8.0;
    fr_b1.dy = (fr_b1.y_max-fr_b1.y_min);
    fit_sub_domain(fr_b1, cheby_order, count_cheby_b1, cheby_coeff_filter_b1,
        cheby_order_0_b1, cheby_order_1_b1, sub_domain_name_b1);
    
    // Generate Chebyshev polynomials for region c (8.0<Y<=20.0)
    std::string sub_domain_name_c = "Sub-Domain C";
    FitRegion fr_c;
    fr_c.num_points = num_points;
    fr_c.side = 1;
    fr_c.y_min = 8.00001;
    fr_c.y_max = 20.0;
    fr_c.dy = (fr_c.y_max-fr_c.y_min);
    fit_sub_domain(fr_c, cheby_order, count_cheby_c, cheby_coeff_filter_c,
        cheby_order_0_c, cheby_order_1_c, sub_domain_name_c);

    // Deallocate heap memory
    mkl_free(cheby_coeff_filter_a0);
    mkl_free(cheby_coeff_filter_a1);
    mkl_free(cheby_coeff_filter_b0);
    mkl_free(cheby_coeff_filter_b1);
    mkl_free(cheby_coeff_filter_c);
    mkl_free(cheby_order_0_a0);
    mkl_free(cheby_order_0_a1);
    mkl_free(cheby_order_0_b0);
    mkl_free(cheby_order_0_b1);
    mkl_free(cheby_order_0_c);
    mkl_free(cheby_order_1_a0);
    mkl_free(cheby_order_1_a1);
    mkl_free(cheby_order_1_b0);
    mkl_free(cheby_order_1_b1);
    mkl_free(cheby_order_1_c);

    return 0;
}


// void write_domain(std::infile infile, )