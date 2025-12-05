
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

// Include local modules
#include "../config.hpp"


class P2
{
public:
    cusfloat*       c;
    int             current_inter = -1;
    int             current_inter_x = -1;
    int             current_inter_y = -1;
    cusfloat*       intervals_bounds_x;
    cusfloat*       intervals_bounds_y;
    bool            is_build = false;
    int*            nx;
    int*            ny;
    int             num_c = 0;
    int             num_intervals_x = 0;
    int             num_intervals_y = 0;
    int*            num_points;
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

    void calculate_interval(cusfloat x, cusfloat y);
    cusfloat calculate_xy_cheby(cusfloat x, cusfloat y);
    cusfloat fit_boundary_x(cusfloat x);
    cusfloat fit_boundary_y(cusfloat x);
    void initialize(void);
    cusfloat x_map(cusfloat x);
    cusfloat y_map(cusfloat y);

};


class R11: public P2
{
public:
    
    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R11A_dX: public P2
{
public:
    
    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R11B_dX: public P2
{
public:
    
    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R12: public P2
{
public:
    
    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R12_dX: public P2
{
public:
    
    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R21: public P2
{
public:

    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R21_dX: public P2
{
public:

    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R22: public P2
{
public:

    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


class R22_dX: public P2
{
public:

    cusfloat calculate_xy(cusfloat x, cusfloat y);
};


// Declare module functions
void set_data_r11(R11* r11);
void set_data_r11a_dx(R11A_dX* r11);
void set_data_r11b_dx(R11B_dX* r11);
void set_data_r12(R12* r11);
void set_data_r12_dx(R12_dX* r11);
void set_data_r21(R21* r11);
void set_data_r21_dx(R21_dX* r11);
void set_data_r22(R22* r11);
void set_data_r22_dx(R22_dX* r11);
