
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
#include "../../containers/mpi_config.hpp"
#include "../../containers/rad_diff_data.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../waves/wave_dispersion_so.hpp"


using RDDQTF = RadDiffData<cuscomplex, RDDQTFConfig>;


void        calculate_qtf_indirect_body_term(
                                                Input*          input,
                                                std::size_t     freq_pos_i,
                                                std::size_t     freq_pos_j,
                                                QTFTypeE        qtf_type,
                                                cuscomplex*     raos_hist,
                                                RDDQTF*         body_gp,
                                                RDDQTF*         wl_gp,
                                                cuscomplex*     qtf_body_force
                                            );


void        calculate_qtf_indirect_fs_near_term(
                                                    Input*          input,
                                                    std::size_t     freq_pos_i,
                                                    std::size_t     freq_pos_j,
                                                    QTFTypeE        qtf_type,
                                                    RDDQTF*         body_gp,
                                                    cuscomplex*     qtf_fs_force
                                                );


void        calculate_qtf_indirect_fs_far_term(
                                                    Input*          input,
                                                    std::size_t     freq_pos_i,
                                                    std::size_t     freq_pos_j,
                                                    QTFTypeE        qtf_type,
                                                    SimulationData* sim_data,
                                                    cuscomplex*     qtf_fs_force
                                                );


cuscomplex  calculate_r0_integral(
                                                    cusfloat            R,
                                                    WaveDispersionSO*   wdso,
                                                    int                 l_order,
                                                    QTFTypeE            qtf_type
                                );


cuscomplex  calculate_r1_integral(
                                                    cusfloat            R,
                                                    WaveDispersionSO*   wdso,
                                                    int                 l_order,
                                                    QTFTypeE            qtf_type
                                );


void  		calculate_theta_integral(
                                                    Input*      input,
                                                    cusfloat    beta,
                                                    int         l_order,
                                                    QTFTypeE    qtf_type,
                                                    int         theta_type,
                                                    cuscomplex* kochin_cos_pert_j,
                                                    cuscomplex* kochin_sin_pert_j,
                                                    cuscomplex* kochin_cos_rad_i,
                                                    cuscomplex* kochin_sin_rad_i,
                                                    cuscomplex* kochin_cos_rad_j,
                                                    cuscomplex* kochin_sin_rad_j,
                                                    cuscomplex* body_force
                                    );


void        calculate_secord_force_indirect(
                                                    Input*          input,
                            MpiConfig*      mpi_config,
                                                    MeshGroup*      mesh_gp,
                            std::size_t     freq_pos_i,
                            std::size_t     freq_pos_j,
                            QTFTypeE        qtf_type,
                            RDDQTF*         body_gp,
                            RDDQTF*         fs_gp,
                            RDDQTF*         wl_gp,
                                                    SimulationData* sim_data
                                            );
