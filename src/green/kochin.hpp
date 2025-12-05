
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
#include "../containers/simulation_data.hpp"
#include "../inout/input.hpp"
#include "../interfaces/kochin_interface.hpp"
#include "../mesh/mesh_group.hpp"


void        calculate_kochin_coefficients(
                                            Input*              input,
                                            MeshGroup*          mesh_gp,
                                            KochinInterface*    kochin,
                                            cuscomplex*         sources,
                                            cuscomplex*         cos_coeff,
                                            cuscomplex*         sin_coeff
                                        );


cusfloat    calculate_kochin_cosexp_t0(
                                            cusfloat            beta,
                                            int                 ln,
                                            int                 m,
                                            int                 n
                                        );


cusfloat    calculate_kochin_cosexp_t1(
                                            cusfloat            beta,
                                            cusfloat            ln,
                                            cusfloat            m,
                                            cusfloat            n
                                        );


cusfloat    calculate_kochin_cosexp_t2(
                                            cusfloat            beta,
                                            cusfloat            ln,
                                            cusfloat            m,
                                            cusfloat            n
                                        );


cusfloat    calculate_kochin_cosexp_t3(
                                            cusfloat            beta,
                                            cusfloat            ln,
                                            cusfloat            m,
                                            cusfloat            n
                                        );


void        calculate_kochin_pert_coeffs(
                                            Input*              input,
                                            MeshGroup*          mesh_gp,
                                            int                 freq_pos,
                                            SimulationData*     sim_data
                                        );


void        calculate_kochin_rad_coeffs(
                                            Input*              input,
                                            MeshGroup*          mesh_gp,
                                            int                 freq_pos,
                                            SimulationData*     sim_data
                                        );
