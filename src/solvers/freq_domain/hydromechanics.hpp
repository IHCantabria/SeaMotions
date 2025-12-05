
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
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
// #include "../../interfaces/hmf_interface.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_hydromechanic_coeffs_lin( 
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panels_pot,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad,
                                                cuscomplex*     panels_pressure
                                            );


// void    calculate_hydromechanic_coeffs_nlin(
//                                                 Input*          input,
//                                                 MpiConfig*      mpi_config,
//                                                 MeshGroup*      mesh_gp,
//                                                 HMFInterface*   hmf_interf,
//                                                 cusfloat        ang_freq,
//                                                 cusfloat*       added_mass,
//                                                 cusfloat*       damping_rad
//                                             );
