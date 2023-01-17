
// Include general usage libraries
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../math/chebyshev.hpp"
#include "../config.hpp"
#include "../math/math_tools.hpp"
#include "pulsating.hpp"
#include "../math/special_math.hpp"

// Redefine-include name spaces
namespace fs = std::filesystem;


// Declare local module functions
struct FitRegion;
void generate_grid(FitRegion fr, double* x, double* y);
void generate_grid_cheby(FitRegion fr, double* x, double* y);
void get_x_range(double y, double &x0, double &x1, int side);
void write_domain(std::ofstream &outfile, int num_cheby, double* cheby_coeff,
                int* cheby_order_0, int* cheby_order_1, std::string domain_key);
template<typename T> void write_vector(std::ofstream &outfile, int num_points, T* vec, int shift);


struct FitRegion
{
    double dy = 0.0;
    int num_points = 0;
    int side = 0;
    double y_min = 0.0;
    double y_max = 0.0;

    void fit(double x, double y, double &xl, double &yl)
    {
        // Calculate x local range
        double x0 = 0.0, x1 = 0.0;
        get_x_range(y, x0, x1, this->side);
        double dx = x1-x0;
        xl = 2.0*(x-x0)/dx-1.0;

        // Calculate y local range
        yl = 2.0*(y-this->y_min)/this->dy-1.0;

    }

};


double eval_chebyshev_xy(int n_order, FitRegion fr, double x, double y, double cheby_prec, double* cheby_coeff)
{
    double sol = 0.0;
    double xl = 0.0, yl = 0.0;
    for (int i=0; i<n_order; i++)
    {
        for (int j=0; j<n_order; j++)
        {
            if (std::abs(cheby_coeff[i*n_order+j])>cheby_prec)
            {
                fr.fit(x, y, xl, yl);
                sol += cheby_coeff[i*n_order+j]*chebyshev_poly_raw(i, xl)*chebyshev_poly_raw(j, yl);
            }
        }
    }

    return sol;
}


void fit_chebyshev_xy(int n_order, int num_points, FitRegion fr, double* x, double * y,
    double* f, double* cheby_coeff)
{
    // Generate kernel matrix
    int num_cols = n_order*n_order;
    double* A = generate_empty_vector<double>(num_points*num_cols);
    double* At = generate_empty_vector<double>(num_points*num_cols);
    double xl = 0.0, yl = 0.0;
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

    double* Ak = generate_empty_vector<double>(num_cols*num_cols);
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


void fit_sub_domain(FitRegion fr, int cheby_order, double cheby_thres, int &count_coeffs,
    double* cheby_coeff_filter, int* cheby_order_0, 
    int* cheby_order_1, std::string sub_domain_name)
{
    // Define allocation matrixes
    int cheby_order_2 = cheby_order*cheby_order;
    int num_points_p2 = fr.num_points*fr.num_points;
    double* x = generate_empty_vector<double>(num_points_p2);
    double* xs = generate_empty_vector<double>(num_points_p2);
    double* y = generate_empty_vector<double>(num_points_p2);
    double* ys = generate_empty_vector<double>(num_points_p2);

    // Generate grid
    generate_grid(fr, xs, ys);
    generate_grid_cheby(fr, x, y);

    // Generate target vector
    double* f = generate_empty_vector<double>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        f[i] = (
            expint_inf_depth_num(x[i], y[i])
            -1/std::sqrt(pow2s(x[i])+pow2s(y[i]))
            +std::exp(-y[i])/x[i]
            )*pow3s(std::sqrt(pow2s(x[i])+pow2s(y[i])))/y[i];
    }

    double* fxy = generate_empty_vector<double>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        fxy[i] = wave_term_inf_depth_num(xs[i], ys[i]) + 1/std::sqrt(pow2s(xs[i])+pow2s(ys[i]));
    }

    // Fit Chebyshev polynomials
    double* cheby_coeff = generate_empty_vector<double>(cheby_order_2);
    fit_chebyshev_xy(cheby_order, num_points_p2, fr, x, y, f, cheby_coeff);

    // Evaluate Chebyshev polynomials
    double* sol_fit = generate_empty_vector<double>(num_points_p2);
    double f1 = 0.0;
    double R = 0.0;
    for (int i=0; i<=num_points_p2; i++)
    {
        R = std::sqrt(pow2s(xs[i])+pow2s(ys[i]));
        f1 = 1/R - std::exp(-ys[i])/xs[i] + ys[i]/pow3s(R)*eval_chebyshev_xy(cheby_order, fr, xs[i], ys[i],
                                                                            cheby_thres, cheby_coeff);
        f1 = 1/R - PI*std::exp(-ys[i])*(struve0(xs[i])+bessely0(xs[i])) - 2.0*f1;
        sol_fit[i] = f1;
    }

    // Compute error with respecto the reference values
    double* err = generate_empty_vector<double>(num_points_p2);
    for (int i=0; i<num_points_p2; i++)
    {
        err[i] = sol_fit[i]-fxy[i];
    }

    // Compute error statistics
    int count_thr = 0;
    double max_err = 0.0;
    double mean_err = 0.0;
    double min_err = 0.0;
    double threshold = 1e-6;
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
    for (int i=0; i<cheby_order_2; i++)
    {
        if (std::abs(cheby_coeff[i])>=cheby_thres)
        {
            count_coeffs++;
        }
    }

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


void generate_grid(FitRegion fr, double* x, double* y)
{
    // Generate y coordinate grid
    double dy = (fr.y_max-fr.y_min)/(fr.num_points-1);
    double yi = 0.0;
    for (int i=0; i<fr.num_points; i++)
    {
        yi = fr.y_min + i*dy;
        for (int j=0; j<fr.num_points; j++)
        {
            y[i*fr.num_points+j] = yi;
        }
    }

    // Generate x coordinate grid
    double dx = 0.0;
    double x0 = 0.0, x1 = 0.0;

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


void generate_grid_cheby(FitRegion fr, double* x, double* y)
{
    // Get Chebyshev polynomial roots
    cusfloat* cheby_roots = generate_empty_vector<cusfloat>(fr.num_points);
    chebyshev_poly_roots(fr.num_points, cheby_roots);

    // Generate y coordinate grid
    double dy = fr.y_max-fr.y_min;
    double yi = 0.0;
    for (int i=0; i<fr.num_points; i++)
    {
        yi = fr.y_min + (cheby_roots[i]+1)*dy/2.0;
        for (int j=0; j<fr.num_points; j++)
        {
            y[i*fr.num_points+j] = yi;
        }
    }

    // Generate x coordinate grid
    double dx = 0.0;
    double x0 = 0.0, x1 = 0.0;

    for (int i=0; i<fr.num_points; i++)
    {
        get_x_range(y[i*fr.num_points], x0, x1, fr.side);
        dx = x1-x0;
        for (int j=0; j<fr.num_points; j++)
        {
            x[i*fr.num_points+j] = x0 + (cheby_roots[j]+1)*dx/2.0;
        }
    }

    // Deallocate heap memory
    mkl_free(cheby_roots);

}


void get_x_range(double y, double &x0, double &x1, int side)
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


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    fs::path folder_path(argv[1]);

    // Define fit properties
    int num_points = 100;
    int cheby_order = 20;
    int cheby_order_2 = cheby_order*cheby_order;

    // Allocate heap memory to storage fit results
    double* cheby_coeff_filter_a0 = generate_empty_vector<double>(cheby_order_2);
    double* cheby_coeff_filter_a1 = generate_empty_vector<double>(cheby_order_2);
    double* cheby_coeff_filter_b0 = generate_empty_vector<double>(cheby_order_2);
    double* cheby_coeff_filter_b1 = generate_empty_vector<double>(cheby_order_2);
    double* cheby_coeff_filter_c = generate_empty_vector<double>(cheby_order_2);
    int* cheby_order_0_a0 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_0_a1 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_0_b0 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_0_b1 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_0_c = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_1_a0 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_1_a1 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_1_b0 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_1_b1 = generate_empty_vector<int>(cheby_order_2);
    int* cheby_order_1_c = generate_empty_vector<int>(cheby_order_2);
    int count_cheby_a0 = 0;
    int count_cheby_a1 = 0;
    int count_cheby_b0 = 0;
    int count_cheby_b1 = 0;
    int count_cheby_c = 0;

    // Generate Chebyshev polynomials for region a0 (X<3.0 & 2.0<=Y<=4.0)
    std::string sub_domain_name_a0 = "Sub-Domain A0";
    double cheby_thres_a0 = 1e-6;
    FitRegion fr_a0;
    fr_a0.num_points = num_points;
    fr_a0.side = 0;
    fr_a0.y_min = 2.0;
    fr_a0.y_max = 4.0;
    fr_a0.dy = (fr_a0.y_max-fr_a0.y_min);
    fit_sub_domain(fr_a0, cheby_order, cheby_thres_a0, count_cheby_a0, cheby_coeff_filter_a0,
        cheby_order_0_a0, cheby_order_1_a0, sub_domain_name_a0);

    // Generate Chebyshev polynomials for region a1 (X>3.0 & 2.0<=Y<=4.0)
    std::string sub_domain_name_a1 = "Sub-Domain A1";
    double cheby_thres_a1 = 1e-5;
    FitRegion fr_a1;
    fr_a1.num_points = num_points;
    fr_a1.side = 1;
    fr_a1.y_min = 2.0;
    fr_a1.y_max = 4.0;
    fr_a1.dy = (fr_a1.y_max-fr_a1.y_min);
    fit_sub_domain(fr_a1, cheby_order, cheby_thres_a1, count_cheby_a1, cheby_coeff_filter_a1,
        cheby_order_0_a1, cheby_order_1_a1, sub_domain_name_a1);

    // Generate Chebyshev polynomials for region b0 (X<8.0 & 4.0<Y<=8.0)
    std::string sub_domain_name_b0 = "Sub-Domain B0";
    double cheby_thres_b0 = 1e-5;
    FitRegion fr_b0;
    fr_b0.num_points = num_points;
    fr_b0.side = 0;
    fr_b0.y_min = 4.0001;
    fr_b0.y_max = 8.0;
    fr_b0.dy = (fr_b0.y_max-fr_b0.y_min);
    fit_sub_domain(fr_b0, cheby_order, cheby_thres_b0, count_cheby_b0, cheby_coeff_filter_b0,
        cheby_order_0_b0, cheby_order_1_b0, sub_domain_name_b0);

    // Generate Chebyshev polynomials for region b1 (X>8.0 & 4.0<Y<=8.0)
    std::string sub_domain_name_b1 = "Sub-Domain B1";
    double cheby_thres_b1 = 1e-5;
    FitRegion fr_b1;
    fr_b1.num_points = num_points;
    fr_b1.side = 1;
    fr_b1.y_min = 4.0001;
    fr_b1.y_max = 8.0;
    fr_b1.dy = (fr_b1.y_max-fr_b1.y_min);
    fit_sub_domain(fr_b1, cheby_order, cheby_thres_b1, count_cheby_b1, cheby_coeff_filter_b1,
        cheby_order_0_b1, cheby_order_1_b1, sub_domain_name_b1);
    
    // Generate Chebyshev polynomials for region c (8.0<Y<=20.0)
    std::string sub_domain_name_c = "Sub-Domain C";
    double cheby_thres_c = 1e-5;
    FitRegion fr_c;
    fr_c.num_points = num_points;
    fr_c.side = 1;
    fr_c.y_min = 8.00001;
    fr_c.y_max = 20.0;
    fr_c.dy = (fr_c.y_max-fr_c.y_min);
    fit_sub_domain(fr_c, cheby_order, cheby_thres_c, count_cheby_c, cheby_coeff_filter_c,
        cheby_order_0_c, cheby_order_1_c, sub_domain_name_c);

    // Write module
    fs::path file_path = folder_path / "chebyshev_inf_depth.hpp";
    std::ofstream outfile(file_path);

    outfile << std::endl;
    outfile << "#ifndef __chebyshev_inf_depth_hpp" << std::endl;
    outfile << "#define __chebyshev_inf_depth_hpp" << std::endl;
    outfile << std::endl << std::endl;

    outfile << "#include \"../config.hpp\"" << std::endl << std::endl;

    outfile << "namespace chebyinf{" << std::endl;

    write_domain(outfile, count_cheby_a0, cheby_coeff_filter_a0, cheby_order_0_a0,
                cheby_order_1_a0, "a0");

    write_domain(outfile, count_cheby_a1, cheby_coeff_filter_a1, cheby_order_0_a1,
                cheby_order_1_a1, "a1");

    write_domain(outfile, count_cheby_b0, cheby_coeff_filter_b0, cheby_order_0_b0,
                cheby_order_1_b0, "b0");

    write_domain(outfile, count_cheby_b1, cheby_coeff_filter_b1, cheby_order_0_b1,
                cheby_order_1_b1, "b1");

    write_domain(outfile, count_cheby_c, cheby_coeff_filter_c, cheby_order_0_c,
                cheby_order_1_c, "c");
    
    outfile << "}" << std::endl;

    outfile << std::endl;
    outfile << "#endif" << std::endl;


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


// void write_domain(std::ofstream &outfile, int num_cheby, double* cheby_coeff_filter, 
//                 double* cheby_order_0, double* cheby_order_1, std::string domain_key)
void write_domain(std::ofstream &outfile, int num_cheby, double* cheby_coeff, int* cheby_order_0,
    int* cheby_order_1, std::string domain_key)
{
    int num_pad_space = 32;
    std::string pad_space(num_pad_space, ' ');
    outfile << "    ///////////////////////////////////////////" << std::endl;
    outfile << "    ////////////// sub_domain " << domain_key <<"///////////////" << std::endl;
    outfile << "    ///////////////////////////////////////////" << std::endl;
    outfile << "    constexpr int num_cheby_" << domain_key << " = " << num_cheby << ";" << std::endl;
    outfile << "    constexpr cusfloat cheby_coeff_" << domain_key << "[" << num_cheby << "] = " << "{" << std::endl;
    write_vector(outfile, num_cheby, cheby_coeff, num_pad_space);
    outfile << pad_space << "};" << std::endl;
    outfile << "    constexpr int cheby_order_0_" << domain_key << "[" << num_cheby << "] = " << "{" << std::endl;
    write_vector(outfile, num_cheby, cheby_order_0, num_pad_space);
    outfile << pad_space << "};" << std::endl;
    outfile << "    constexpr int cheby_order_1_" << domain_key  << "[" << num_cheby << "] = " << "{" << std::endl;
    write_vector(outfile, num_cheby, cheby_order_1, num_pad_space);
    outfile << pad_space << "};" << std::endl;

    outfile << std::endl << std::endl;
}

template<typename T>
void write_vector(std::ofstream &outfile, int num_points, T* vec, int shift)
{
    std::string shift_space(shift, ' ');
    for (int i=0; i<num_points; i++)
    {
        outfile << shift_space;
        outfile << std::setprecision(16) << std::scientific << vec[i];
        outfile << ", // C[" << i << "]" << std::endl;
    }
}