
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
#include "../config.hpp"
#include "../containers/containers.hpp"


void    calculate_source_monopole_potential(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat&       phi
                                            );


void    calculate_source_monopole_velocity(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat*       field_point,
                                            cusfloat*       velocity
                                        );


void    calculate_source_potential_newman(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag, 
                                            int             multipole_flag, 
                                            cusfloat&       phi
                                        );


void    calculate_source_velocity_newman(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag,
                                            int             multipole_flag, 
                                            cusfloat        *velocity
                                        );


void    calculate_source_newman(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag,
                                            int             multipole_flag, 
                                            cusfloat        *velocity,
                                            cusfloat&       phi
                                        );


void    calculate_source_potential_hess(
                                            PanelGeom*      panel,
                                            cusfloat*       field_point,
                                            int             fp_local_flag,
                                            int             multipole_flag,
                                            cusfloat        &phi
                                        );


void    calculate_source_velocity_hess(
                                            PanelGeom*      panel,
                                            cusfloat*       field_point, 
                                            int             fp_local_flag, 
                                            int             multipole_flag,
                                            cusfloat        *velocity
                                        );
