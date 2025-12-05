
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


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