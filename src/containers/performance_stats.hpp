
#ifndef __performance_tests_hpp
#define __performance_tests_hpp

struct PerformanceStats
{
    int count = 0;
    double t_max = 0.0;
    double t_mean = 0.0;
    double t_min = 1e16;

    void add_performance_point( double t );

};

#endif