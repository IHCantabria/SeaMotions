
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
#include "../containers/panel_geom.hpp"
#include "../math/custensor/custensor.hpp"
#include "rad_diff_data_config.hpp"


template<typename T, typename Config>
struct PanelData
{
    /*** Template parameters from Config ***/
    static constexpr int mode_comp      = Config::mode_comp;
    static constexpr int mode_f         = Config::mode_f;
    static constexpr int mode_dfdn      = Config::mode_dfdn;
    static constexpr int mode_dfdc      = Config::mode_dfdc;
    static constexpr int store_modes    = Config::store_modes;
    static constexpr int store_freqs    = Config::store_freqs;

private:
    // Define private variables
    bool                            _is_heap        = false;    // Flag to indicate if memory is allocated on heap

    /* Declare class private methods */
    void    _allocate_memory( 
                                    std::size_t field_points_np_,
                                    std::size_t freqs_np_,
                                    std::size_t headings_np_,
                                    std::size_t dofs_np_
                                );

    void    _load_field_points( 
                                    PanelGeom*  panel_geom_,
                                    bool        use_waterline_
                                );

public:
    // Declare public variables
    std::size_t                     body_id         = 0;        // Index of the body at which the panel belongs
    std::size_t                     dofs_np         = 0;        // Number of degrees of freedom
    std::size_t                     field_points_np = 0;        // Number of field points
    cut::CusTensor<cusfloat>        field_points;               // Store field points coordinates
    std::size_t                     freqs_np        = 0;        // Number of frequencies
    std::size_t                     headings_np     = 0;        // Number of headings
    PanelGeom*                      panel_geom      ;           // Pointer to panel geometry information
    cut::CusTensor<T>               pot_incident    ;           // Store wave incident potential value                           [frequencies, headings, field_points]
    cut::CusTensor<T>               pot_rad         ;           // Store wave radiated potential value                           [frequencies, headings, field_points]
    cut::CusTensor<T>               pot_diff        ;           // Store wave diffracted potential value                         [frequencies, headings, field_points]
    cut::CusTensor<T>               pot_total       ;           // Store wave total potential value                              [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               pot_total_2     ;           // Store wave total potential value                              [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               pot_modes       ;           // Store raw radiation-diffraction potential modes               [frequencies, dofs + headings, field_points]
    cut::CusTensor<T>               press_incident  ;           // Store incident pressure value                                    [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               press_rad       ;           // Store radiated pressure value                                    [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               press_diff      ;           // Store diffracted pressure value                                    [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               press_total     ;           // Store total pressure value                                    [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               vel_dn_incident ;           // Store wave incident normal velocity derivative value          [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_dn_rad      ;           // Store radiated normal velocity derivative values              [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_dn_diff     ;           // Store diffracted normal velocity derivative values            [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_dn_total    ;           // Store total normal velocity derivative value                  [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               vel_x_incident  ;           // Store wave incident velocity_x value                          [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_y_incident  ;           // Store wave incident velocity_y value                          [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_z_incident  ;           // Store wave incident velocity_z value                          [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_x_rad       ;           // Store radiated velocity_x values                              [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_y_rad       ;           // Store radiated velocity_y values                              [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_z_rad       ;           // Store radiated velocity_z values                              [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_x_diff      ;           // Store diffracted velocity_x values                            [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_y_diff      ;           // Store diffracted velocity_y values                            [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_z_diff      ;           // Store diffracted velocity_z values                            [frequencies, headings, field_points]
    cut::CusTensor<T>               vel_x_total     ;           // Store total velocity_x value                                  [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               vel_y_total     ;           // Store total velocity_y value                                  [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               vel_z_total     ;           // Store total velocity_z value                                  [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               vel_x_modes     ;           // Store raw radiation-diffraction velocity_x modes              [frequencies, dofs + headings, field_points]
    cut::CusTensor<T>               vel_y_modes     ;           // Store raw radiation-diffraction velocity_y modes              [frequencies, dofs + headings, field_points]
    cut::CusTensor<T>               vel_z_modes     ;           // Store raw radiation-diffraction velocity_z modes              [frequencies, dofs + headings, field_points]
    cut::CusTensor<T>               wev_total       ;           // Store wave elevation                                          [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }
    cut::CusTensor<T>               wev_rel_total   ;           // Store relative wave elevation                                 [frequencies, headings, field_points]             | Total composition: { incident + sum( epsi · radiation ) + diffraction }


    /* Declare class contructor */
    PanelData( ) = default;

    PanelData( 
                std::size_t field_points_np_,
                cusfloat*   field_points_,
                std::size_t freqs_np_,
                std::size_t headings_np_,
                std::size_t dofs_np_
            );

    PanelData( 
                PanelGeom*  panel_geom_,
                std::size_t body_id_,
                std::size_t freqs_np_,
                std::size_t headings_np_,
                std::size_t dofs_np_,
                bool        use_waterline_ = false
            );

    /* Declare public methods */
    void clear_data( void );

    template<FieldTypeE field_type, FieldComponentE field_component>
    const cut::CusTensor<T>* get_field_data( void ) const;

};


// Include template implementation file
#include "panel_data.txx"