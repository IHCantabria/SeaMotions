
// Include local modules
#include "performance_stats.hpp"


void PerformanceStats::add_performance_point(double t)
{
    // Increment counter to have this new iteration
    this->count++;

    // Check if the new performance point is higher than the maximum
    if (t > t_max)
    {
        this->t_max = t;
    }

    // Check if the new performance point is lower than the minimum
    if (t < t_min)
    {
        this->t_min = t;
    }

    // Add performance point to the mean
    this->t_mean = (this->t_mean*(this->count-1)+t)/this->count;
}