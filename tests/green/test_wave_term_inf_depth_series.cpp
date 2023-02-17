
// Include general usage libraries
#include <functional>
#include <iomanip>
#include <iostream>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating_inf_depth_series.hpp"
#include "../../src/green/pulsating_inf_depth_cheby.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/special_math.hpp"


void calculate_error_stats(int N, cusfloat* err, cusfloat threshold, int &count_threshold,
    cusfloat &max_err, cusfloat &mean_err, cusfloat &min_err)
{
    // Initialize statistic values
    count_threshold = 0;
    max_err = 0.0;
    mean_err = 0.0;
    min_err = 1e16;

    // Loop over error vector to calculate statistics
    for (int i=0; i<N; i++)
    {
        // Check if the value goes over threshold
        if (err[i] > threshold)
        {
            count_threshold++;
        }

        // Check maximum value
        if (err[i] > max_err)
        {
            max_err = err[i];
        }

        // Check minimum value
        if (err[i] < min_err)
        {
            min_err = err[i];
        }

        // Add value to account for the mean 
        mean_err += err[i];
    }
    mean_err /= N;

}


void generate_domain(int N, cusfloat* X, cusfloat* Y, cusfloat x_min, cusfloat x_max,
    cusfloat y_min, cusfloat y_max)
{
    // Generate X coordinate points
    cusfloat dx = (x_max-x_min)/(N-1);
    for (int i=0; i<N; i++)
    {
        X[i] = x_min + i*dx;
    }

    // Generate Y coordinate points
    cusfloat dy = (y_max-y_min)/(N-1);
    for (int i=0; i<N; i++)
    {
        Y[i] = y_min + i*dy;
    }
}


bool launch_test(int N, cusfloat* X, cusfloat* Y, std::function<cusfloat(cusfloat,cusfloat)> f_num,
    std::function<cusfloat(cusfloat,cusfloat)> f_ser, bool show_stats_force)
{
    // Allocate space for error statistics
    cusfloat* err = generate_empty_vector<cusfloat>(N*N);

    // Loop over grid points to check the accuracy of the model
    bool pass = true;
    cusfloat diff = 0.0, diff_log = 0.0;
    cusfloat f_ref = 0.0, f_ref_log = 0.0, f_mod = 0.0;
    int count = 0;
    cusfloat threshold = 3e-6;
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Calculate reference value using Romberg integration method
            f_ref = f_num(X[i], Y[j]);
            f_ref_log = std::log10(std::abs(f_ref));

            // Calculate model value
            f_mod = f_ser(X[i], Y[j]);

            // Evaluate difference
            diff = f_mod - f_ref;
            diff_log = std::log10(std::abs(diff));
            if (((f_ref_log - diff_log)<6) && (std::abs(diff)>threshold))
            {
                err[count] = std::abs(diff);
                count++;
                // std::cerr << std::setprecision(6) << "X: " << X[i] << " - Y: " << Y[j] << " - f_ref: " << f_ref;
                // std::cerr << " - f_mod: " << f_mod << " - diff: " << diff << std::endl;
                // std::cout << "--> f_ref_log: " << f_ref_log << " - diff_log: " << diff_log << " - diff_order: " << (f_ref_log - diff_log) << std::endl;
            }
        }
    }

    // Calculate error statistics
    int count_threshold = 0;
    cusfloat max_err = 0.0, mean_err = 0.0, min_err = 0.0;
    calculate_error_stats(count, err, threshold, count_threshold, max_err, mean_err, min_err);

    if ((mean_err>7.0e-6) || (max_err>3.5e-5) || show_stats_force)
    {
        std::cout << "Test error statistics:" << std::endl;
        std::cout << " - Max. Error: " << max_err << std::endl;
        std::cout << " - Min. Error: " << min_err << std::endl;
        std::cout << " - Mean Error: " << mean_err << std::endl;
        std::cout << " - Test threshold: " << threshold << " - Count Over threshold: " << count_threshold << std::endl;
    }

    // Deallocate heap memory
    mkl_free(err);

    return pass;
}


int main(void)
{
    // Declare local variables
    int pass = false;

    // Generate computational grid
    constexpr int N = 1000;
    cusfloat X[N], Y[N];
    generate_domain(N, X, Y, 0.001, 40.0, 0.001, 40.0);

    // Test modelling function over all XY plane
    pass = launch_test(
        N, 
        X, 
        Y, 
        wave_term_inf_depth_num,
        wave_term_inf_depth_series,
        false
        );
    if (!pass)
    {
        std::cerr << "test_wave_term_inf_depth_series failed!" << std::endl;
        return 1;
    }

    // Test modelling function horizontal derivative over all XY plane
    pass = launch_test(
        N, 
        X, 
        Y, 
        wave_term_inf_depth_num_dxndim,
        wave_term_inf_depth_dxndim_series,
        false
        );
    if (!pass)
    {
        std::cerr << "test_wave_term_inf_depth_dxndim_series failed!" << std::endl;
        return 1;
    }

    // Test modelling function vertical derivative over all XY plane
    pass = launch_test(
        N, 
        X, 
        Y, 
        wave_term_inf_depth_num_dyndim,
        wave_term_inf_depth_dyndim_series,
        false
        );
    if (!pass)
    {
        std::cerr << "test_wave_term_inf_depth_dyndim_series failed!" << std::endl;
        return 1;
    }

    return 0;
}