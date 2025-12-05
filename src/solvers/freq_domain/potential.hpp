
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
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_influence_potmat_steady(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        pot_gp
                                            );


void    calculate_influence_potmat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        pot_gp
                                    );


void    calculate_potpanel_raddif_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     intensities,
                                                cusfloat        ang_freq
                                        );


void    calculate_potpanel_total_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     intensities,
                                                cuscomplex*     raos,
                                                MLGCmpx*        pot_gp,
                                                cuscomplex*     potpanel_total
                                    );
