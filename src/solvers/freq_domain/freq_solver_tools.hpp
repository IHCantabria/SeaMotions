
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
// #include "../../interfaces/hmf_interface.hpp"
// #include "../../interfaces/grfdn_interface.hpp"
// #include "../../interfaces/gwfdn_interface.hpp"
#include "../../hydrostatics.hpp"
#include "../../inout/input.hpp"
#include "../../inout/output.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_freq_domain_coeffs(
                                                MpiConfig*      mpi_config,
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                    );


void    calculate_global_hydstiffness(
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                cusfloat*       hydstiffness
                                    );


void    calculate_global_structural_mass(
                                                Input*          input,
                                                cusfloat*       structural_mass_p0
                                        );


void    freq_domain_linear_solver(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                MeshGroup*      mesh_fs_qtf_gp,
                                                SclCmpx*        scl,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                );


void    freq_domain_nonlinear_solver(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                SclCmpx*        scl,
                                                Hydrostatics**  hydrostatics,
                                                Output*         output
                                    );
