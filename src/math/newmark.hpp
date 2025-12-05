
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
#include "math_tools.hpp"
#include "./sparse/sparse_math.hpp"


template<typename T>
class NewmarkBeta
{
public:
    // Define class attributes
    cusfloat            beta            = 0.0;
    int                 count           = 0;
    sparse_matrix_t*    damp_mat        = nullptr;
    T                   fext_fcn          ;
    cusfloat            gamma           = 0.0;
    CSRMatrix*          lhs_csr         = nullptr;
    SparseSolver*       linsolve        = nullptr;
    sparse_matrix_t*    mass_mat        = nullptr;
    int                 mat_size        = 0;
    sparse_matrix_t*    stiff_damp_acc  = nullptr;
    sparse_matrix_t*    stiff_damp_vel  = nullptr;
    sparse_matrix_t*    stiff_mat       = nullptr;
    cusfloat*           rhs             = nullptr;      
    int                 rows_np         = 0;
    cusfloat            time            = 0.0;
    cusfloat            time_init       = 0.0;
    cusfloat            time_old        = 0.0;
    cusfloat            time_step       = 0.0;
    cusfloat*           _u_beta         = nullptr;
    cusfloat*           _u_gamma        = nullptr;
    cusfloat*           y0_acc          = nullptr;
    cusfloat*           y0_pos          = nullptr;
    cusfloat*           y0_vel          = nullptr;
    cusfloat*           y_acc           = nullptr;
    cusfloat*           y_pos           = nullptr;
    cusfloat*           y_vel           = nullptr;
    cusfloat*           y_acc_old       = nullptr;
    cusfloat*           y_pos_old       = nullptr;
    cusfloat*           y_vel_old       = nullptr;

    // Define class constructors and destructors
    NewmarkBeta(
                                        T               fext_in,
                                        CSRMatrix*      mass_mat_in,
                                        CSRMatrix*      stiff_mat_in,
                                        CSRMatrix*      damp_mat_in,
                                        cusfloat        time_step_in,
                                        cusfloat        t0_in,
                                        cusfloat*       y0_pos_in,
                                        cusfloat*       y0_vel_in,
                                        cusfloat*       y0_acc_in
                );
    
    ~NewmarkBeta(void);

    // Define class methods
    void    get_values_at(
                                        cusfloat        time,
                                        cusfloat*       y_pos_at,
                                        cusfloat*       y_vel_at,
                                        cusfloat*       y_acc_at
                        );
    
    void    initialize(
                                        void
                        );
    
    void    interpolate_solution(
                                        cusfloat        time, 
                                        cusfloat*       y_pos_itp,
                                        cusfloat*       y_vel_itp,
                                        cusfloat*       y_acc_itp
                                );

    void step(
                                        void
            );

};


#include "newmark.txx"
