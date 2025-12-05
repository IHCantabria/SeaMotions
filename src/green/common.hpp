
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


void    calculate_distance_node_field(
                                        PanelGeom*  panel, 
                                        cusfloat*   field_point_local, 
                                        cusfloat*   node_fieldp_mod, 
                                        cusfloat*   node_fieldp_dx, 
                                        cusfloat*   node_fieldp_dy, 
                                        cusfloat*   node_fieldp_dz
                                    );
void    calculate_nodes_distance(
                                        PanelGeom*  panel, 
                                        cusfloat*   delta_xi, 
                                        cusfloat*   delta_eta
                                );
void    calculate_polar_angles(
                                        PanelGeom*  panel, 
                                        cusfloat*   delta_xi, 
                                        cusfloat*   delta_eta, 
                                        cusfloat*   polar_angles
                                );
void    calculate_sides_len_local(
                                        PanelGeom*  panel, 
                                        cusfloat*   delta_xi, 
                                        cusfloat*   delta_eta, 
                                        cusfloat* sides_len
                                    );
