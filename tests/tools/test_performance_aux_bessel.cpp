
// Include general usage external libraries
#include <iostream>

// Include local modules
#include "../../src/containers/containers.hpp"
#include "../../src/math/special_math.hpp"
#include "../../src/tools.hpp"


// Declare functions in order to have access along the module
double launch_test(cusfloat (*f)(cusfloat), int num_reps, cusfloat a, cusfloat b);
void print_stats(double pf0, double pf1, double pth0, double pth1, double rf0, double rf1, 
                double rth0, double rth1);


double launch_test(cusfloat (*f)(cusfloat), int num_reps, cusfloat a, cusfloat b)
{
    // Define local variables
    cusfloat m = (b-a)/num_reps;
    cusfloat x = 0.0;

    // Perform iterations
    double t0 = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        x = a + m*i;
        f(x);
    }
    double t1 = get_cpu_time();

    return t1-t0;
}


int main(void)
{
    // Define test domain boundaries
    cusfloat a = 0.01;
    cusfloat b = 1000;

    // Define statistical parameters
    int num_reps = 1e8;
    int num_stats_reps = 10;

    // Create storage objects
    PerformanceStats pf0_stats;
    PerformanceStats pf1_stats;
    PerformanceStats pth0_stats;
    PerformanceStats pth1_stats;
    PerformanceStats rf0_stats;
    PerformanceStats rf1_stats;
    PerformanceStats rth0_stats;
    PerformanceStats rth1_stats;

    // Loop over statistic values
    for (int i=0; i<num_stats_reps; i++)
    {
        // Perform test polynomial_f0
        pf0_stats.add_performance_point(launch_test(polynomial_f0, num_reps, a, b));

        // Perform test polynomial_f1
        pf1_stats.add_performance_point(launch_test(polynomial_f1, num_reps, a, b));

        // Perform test polynomial_th0
        pth0_stats.add_performance_point(launch_test(polynomial_th0, num_reps, a, b));

        // Perform test polynomial_th1
        pth1_stats.add_performance_point(launch_test(polynomial_th1, num_reps, a, b));

        // Perform test rational_fraction_f0
        rf0_stats.add_performance_point(launch_test(rational_fraction_f0, num_reps, a, b));

        // Perform test rational_fraction_f1
        rf1_stats.add_performance_point(launch_test(rational_fraction_f1, num_reps, a, b));

        // Perform test rational_fraction_th0
        rth0_stats.add_performance_point(launch_test(rational_fraction_th0, num_reps, a, b));

        // Perform test rational_fraction_th1
        rth1_stats.add_performance_point(launch_test(rational_fraction_th1, num_reps, a, b));

    }

    // Print-out statistics
    std::cout << "MAXIMUM VALUES:" << std::endl;
    print_stats(pf0_stats.t_max, pf1_stats.t_max, pth0_stats.t_max, pth1_stats.t_max,
                rf0_stats.t_max, rf1_stats.t_max, rth0_stats.t_max, rth1_stats.t_max);
    std::cout << std::endl;
    std::cout << "MINIMUM VALUES:" << std::endl;
    print_stats(pf0_stats.t_min, pf1_stats.t_min, pth0_stats.t_min, pth1_stats.t_min,
                rf0_stats.t_min, rf1_stats.t_min, rth0_stats.t_min, rth1_stats.t_min);
    std::cout << std::endl;
    std::cout << "MEAN VALUES:" << std::endl;
    print_stats(pf0_stats.t_mean, pf1_stats.t_mean, pth0_stats.t_mean, pth1_stats.t_mean,
                rf0_stats.t_mean, rf1_stats.t_mean, rth0_stats.t_mean, rth1_stats.t_mean);
    std::cout << std::endl;

    return 0;
}


void print_stats(double pf0, double pf1, double pth0, double pth1, double rf0, double rf1, 
                double rth0, double rth1)
{
    std::cout << "f0 test -> polynomial: " << pf0 << " - rational fraction: " << rf0 << std::endl;
    std::cout << "f1 test -> polynomial: " << pf1 << " - rational fraction: " << rf1 << std::endl;
    std::cout << "pth0 test -> polynomial: " << pth0 << " - rational fraction: " << rth0 << std::endl;
    std::cout << "pth1 test -> polynomial: " << pth1 << " - rational fraction: " << rth1 << std::endl;
}


// void test_poly_f0()