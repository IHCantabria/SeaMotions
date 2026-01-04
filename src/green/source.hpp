
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


/**
 * @brief Calculate the monopole contribution to velocity and potential
 *
 * @tparam      mode_f          Switch to for the function calculation.
 * @tparam      mode_dfdr       Switch to for the radial derivative calculation.
 * @tparam      mode_dfdz       Switch to for the vertical derivative calculation.
 * @param[in]   panel           Pointer to the panel geometry structure.
 * @param[in]   r0              Distance from the panel center to the field point.
 * @param[in]   field_point     Field point coordinates in local panel system.
 * @param[out]  velocity        Velocity vector induced at the field point.
 * @param[out]  phi             Potential induced at the field point.
 */
template<int mode_f, int mode_dfdr, int mode_dfdz>
void    calculate_source_monopole(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat*       field_point,
                                            cusfloat*       velocity,
                                            cusfloat&       phi
                                        );


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


/**
 * @brief Calculate the multipole contribution to velocity and potential.
 * 
 * @tparam      mode_f          Switch to for the function calculation.
 * @tparam      mode_dfdr       Switch to for the radial derivative calculation.
 * @tparam      mode_dfdz       Switch to for the vertical derivative calculation.
 * @param[in]   panel           Pointer to the panel geometry structure.
 * @param[in]   r0              Distance from the panel center to the field point.
 * @param[in]   field_point     Field point coordinates in local panel system.
 * @param[out]  velocity        Velocity vector induced at the field point.
 * @param[out]  phi             Potential induced at the field point.
 */
template<int mode_f, int mode_dfdr, int mode_dfdz>
void    calculate_source_multipole(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat*       field_point,
                                            cusfloat*       velocity,
                                            cusfloat&       phi
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


/**
 * @brief Calculate the velocity and potential induced by a source panel using
 *        Newman's formulation with templated modes for function and derivatives.
 * 
 * @tparam      mode_f          Switch to for the function calculation.
 * @tparam      mode_dfdr       Switch to for the radial derivative calculation.
 * @tparam      mode_dfdz       Switch to for the vertical derivative calculation.
 * @param[in]   panel           Pointer to the panel geometry structure.
 * @param[in]   field_point     Field point coordinates.
 * @param[in]   fp_local_flag   Flag to indicate if the field point is given in
 *                              global (0) or local (1) coordinates.
 * @param[in]   multipole_flag  Flag to indicate if multipole expansion is allowed.
 * @param[out]  velocity        Velocity vector induced at the field point.
 * @param[out]  phi             Potential induced at the field point.
 */
template<int mode_f, int mode_dfdr, int mode_dfdz>
void    calculate_source_newman_t(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag,
                                            int             multipole_flag, 
                                            cusfloat*       velocity,
                                            cusfloat&       phi
                                        );


// Include template definitions
#include "source.txx"
