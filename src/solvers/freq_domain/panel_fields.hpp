
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

// Function declaration
template<typename T, typename V>    inline  void    calculate_fields_lin(
                                                                                    Input*          input,
                                                                                    MpiConfig*      mpi_config,
                                                                                    MeshGroup*      mesh_gp,
                                                                                    T*              field_func,
                                                                                    V*              fk_fcn,
                                                                                    cusfloat        ang_freq,
                                                                                    cuscomplex*     intensities,
                                                                                    cuscomplex*     raos,
                                                                                    MLGCmpx*        field_gp,
                                                                                    cuscomplex*     field_fk_p0,
                                                                                    cuscomplex*     field_raddif_p0,
                                                                                    cuscomplex*     field_total
                                                                        );

                                    inline  void    calculate_fields_raddif_lin(
                                                                                    Input*          input,
                                                                                    cuscomplex*     intensities,
                                                                                    MLGCmpx*        pot_gp
                                                                                );


template<typename T>                inline  void    calculate_influence_field_mat(
                                                                                    Input*          input,
                                                                                    MeshGroup*      mesh_gp,
                                                                                    T*              field_funct,
                                                                                    MLGCmpx*        field_gp
                                                                                );


                                    inline  void    calculate_pertubation_field(
                                                                                    Input*          input,
                                                                                    cuscomplex*     raos,
                                                                                    cusfloat        ang_freq,
                                                                                    MLGCmpx*        pot_gp,
                                                                                    cuscomplex*     pot_raddif,
                                                                                    cuscomplex*     pot_total
                                                                                );


                                    inline  void    calculate_raddiation_field(
                                                                                    Input*          input,
                                                                                    cuscomplex*     raos,
                                                                                    cusfloat        ang_freq,
                                                                                    MLGCmpx*        pot_gp,
                                                                                    cuscomplex*     pot_raddif,
                                                                                    cuscomplex*     pot_total
                                                                                );


                                    inline  void    calculate_total_field(
                                                                                    Input*          input,
                                                                                    cusfloat        ang_freq,
                                                                                    MLGCmpx*        pot_gp,
                                                                                    cuscomplex*     raos,
                                                                                    cuscomplex*     pot_fk,
                                                                                    cuscomplex*     pot_raddif,
                                                                                    cuscomplex*     pot_total
                                                                            );

// Include module function definition
#include "panel_fields.txx"
