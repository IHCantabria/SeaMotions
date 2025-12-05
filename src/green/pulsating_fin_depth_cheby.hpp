
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

#pragma once

// Include general usage libraries
#include <cassert>
#include <iostream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../config.hpp"
#include "../math/chebyshev.hpp"
#include "./fin_depth_coeffs/L1.hpp"
#include "./fin_depth_coeffs/L1_dA.hpp"
#include "./fin_depth_coeffs/L1_dB.hpp"
#include "./fin_depth_coeffs/L2.hpp"
#include "./fin_depth_coeffs/L3.hpp"
#include "./fin_depth_coeffs/L3_dA.hpp"
#include "./fin_depth_coeffs/L3_dB.hpp"
#include "./fin_depth_coeffs/M1.hpp"
#include "./fin_depth_coeffs/M1_dA.hpp"
#include "./fin_depth_coeffs/M1_dB.hpp"
#include "./fin_depth_coeffs/M2.hpp"
#include "./fin_depth_coeffs/M3.hpp"
#include "./fin_depth_coeffs/M3_dA.hpp"
#include "./fin_depth_coeffs/M3_dB.hpp"


class P3
{
public:
    cusfloat*       c;
    cusfloat*       cf;
    cusfloat*       cf2;
    int             current_inter = -1;
    int             dims = -1;
    cusfloat        int_1d = 0.0;
    cusfloat*       intervals_bounds;
    int*            nx;
    int*            nxf;
    int*            nxf2;
    int*            ny;
    int*            nyf;
    int*            nz;
    int             num_intervals = 0;
    int*            num_points;
    int             num_points_f = 0;
    int             num_points_f2 = 0;
    int*            num_points_cum;
    bool*           x_log_scale;
    cusfloat*       x_max;
    cusfloat*       x_map_scale;
    cusfloat*       x_map_scale_log;
    cusfloat*       x_min;
    cusfloat*       x_min_l10;
    bool*           y_log_scale;
    cusfloat*       y_max;
    cusfloat*       y_map_scale;
    cusfloat*       y_map_scale_log;
    cusfloat*       y_min;
    cusfloat*       y_min_l10;
    bool*           z_log_scale;
    cusfloat*       z_max;
    cusfloat*       z_map_scale;
    cusfloat*       z_map_scale_log;
    cusfloat*       z_min;
    cusfloat*       z_min_l10;

    void calculate_h_1D(cusfloat H);
    void initialize(void);
    void fold_b(cusfloat b);
    void fold_h(cusfloat H);
    int get_interval_h(cusfloat H);
    cusfloat get_value_a(cusfloat a);
    cusfloat get_value_ab(cusfloat a, cusfloat b);
    cusfloat get_value_abh(cusfloat a, cusfloat b, cusfloat h);
    cusfloat x_map(cusfloat x);
    cusfloat y_map(cusfloat y);
    cusfloat z_map(cusfloat z);

};


void set_data_l1(P3 *l1);
void set_data_l1_da(P3 *l1_da);
void set_data_l1_db(P3 *l1_db);
void set_data_l2(P3 *l2);
void set_data_l3(P3 *l3);
void set_data_l3_da(P3 *l3_da);
void set_data_l3_db(P3 *l3_db);
void set_data_m1(P3 *m1);
void set_data_m1_da(P3 *m1_da);
void set_data_m1_db(P3 *m1_db);
void set_data_m2(P3 *m2);
void set_data_m3(P3 *m3);
void set_data_m3_da(P3 *m3_da);
void set_data_m3_db(P3 *m3_db);
