
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
#include "../../inout/input.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_froude_krylov_fo(
                                    Input*          input,
                                    MpiConfig*      mpi_config,
                                    MeshGroup*      mesh_gp,
                                    cusfloat        ang_freq,
                                    cuscomplex*     froude_krylov
                                );

void    calculate_froude_krylov_so(
                                    Input*          input,
                                    MeshGroup*      mesh_gp,
                                    cusfloat        ang_freq_i,
                                    cusfloat        ang_freq_j,
                                    int             qtf_type,
                                    cuscomplex*     froude_krylov
                                );
