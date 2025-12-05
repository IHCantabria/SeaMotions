
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
#include "../../interfaces/grfdn_interface.hpp"
#include "../../interfaces/gwfdn_interface.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_gf_intensity_sysmat(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    GWFDnInterface* gwf_interf,
                                                    cusfloat        w,
                                                    cuscomplex*     sysmat_steady,
                                                    cuscomplex*     sysmat,
                                                    cuscomplex*     sources_int
                                        );


void    calculate_gf_intensity_steady_sysmat_lin(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    cuscomplex*     sysmat
                                                );


void    calculate_gf_intensity_steady_sysmat_nlin(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    GRFDnInterface* grf_interf,
                                                    cuscomplex*     sysmat
                                                );
