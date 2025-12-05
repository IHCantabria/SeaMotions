
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
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void        calculate_field_point_rot(
                                                cuscomplex*     raos_trans,
                                                cuscomplex*     raos_rot,
                                                cusfloat*       field_point,
                                                cusfloat*       cog,
                                                cuscomplex*     point_disp
                                    );


void        calculate_field_point_rot_jac(
                                                cuscomplex*     raos_rot,
                                                cuscomplex*     jac
                                        );


void        calculate_field_point_vel_rot(
                                                cuscomplex*     raos_trans,
                                                cuscomplex*     raos_rot,
                                                cusfloat*       field_point,
                                                cusfloat*       cog,
                                                cusfloat        ang_freq,
                                                cuscomplex*     point_disp
                                        );


std::string compose_dof_path( 
                                                std::string     base_path,
                                                int             dofs_num,
                                                int             ang_freq_num
                            );


void        define_gauss_points_diffrac_panels(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        mat_gp
                                            );


void        define_gauss_points_wl(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        mat_gp
                                );


void        storage_radiation_potential( 
                                                std::string     dof_base_path,
                                                int             field_points_np,
                                                cuscomplex*     rad_potential
                                        );
