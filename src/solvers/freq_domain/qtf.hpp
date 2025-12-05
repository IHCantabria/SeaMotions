
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
#include "../../config.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"

// Set module definitions
#define QTF_DIFF_CODE   0
#define QTF_SUM_CODE    1


// Declare module functions
void        calculate_pinkster(
                                            Input*      input,
                                            MpiConfig*  mpi_config,
                                            MeshGroup*  mesh_gp,
                                            cusfloat    ang_freq_i,
                                            cusfloat    ang_freq_j,
                                            cuscomplex* qtf_values
                            );

cuscomplex  calculate_qtf_diff_term(
                                            cuscomplex c0,
                                            cuscomplex c1
                                    );

void        calculate_qtf_terms_force(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            int             qtf_type,
                                            cuscomplex*     mdrift_rel_we_i,
                                            cuscomplex*     mdrift_rel_we_j,
                                            cuscomplex*     raos_i,
                                            cuscomplex*     raos_j,
                                            cuscomplex*     vel_x_i,
                                            cuscomplex*     vel_y_i,
                                            cuscomplex*     vel_z_i,
                                            cuscomplex*     vel_x_j,
                                            cuscomplex*     vel_y_j,
                                            cuscomplex*     vel_z_j,
                                            cusfloat        ang_freq_i,
                                            cusfloat        ang_freq_j,
                                            cuscomplex*     qtf_values,
                                            cuscomplex*     qtf_wl,
                                            cuscomplex*     qtf_bern,
                                            cuscomplex*     qtf_acc,
                                            cuscomplex*     qtf_mom,
                                            MLGCmpx*        pot_gp,
                                            MLGCmpx*        vel_gp,
                                            bool            is_multi_head
                                    );

void        qtf_distribute_matrix_data(
                                            Input*      input,
                                            int         freq_idx,
                                            int         freq_jdx,
                                            cuscomplex* local_mat,
                                            cuscomplex* global_mat,
                                            int         vector_mode,
                                            int         op_mode
                                        );
