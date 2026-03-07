
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
#include "../containers/morison_coeff_curve.hpp"
#include "../containers/morison_element_def.hpp"
#include "../containers/rad_diff_data.hpp"
#include "../inout/input.hpp"
#include "../math/custensor/custensor.hpp"


class MorisonElement
{
private:
    /*** Define private class attributes ***/
    cusfloat                                    _area_axial         = 0.0;      // Area of the axial component of the Morison element
    cusfloat                                    _area_trans         = 0.0;      // Area of the transverse component of the Morison element
    cusfloat                                    _area_vert          = 0.0;      // Area of the
    MorisonCoeffCurve*                          _ca_axial           = nullptr;  // Added mass coefficient for the axial component of the Morison element
    MorisonCoeffCurve*                          _ca_trans           = nullptr;  // Added mass coefficient for the transverse component of the Morison element
    MorisonCoeffCurve*                          _ca_vert            = nullptr;  // Added mass coefficient for the vertical component of the Morison element
    cut::CusTensor<cusfloat>                    _field_points_l     ;           // Field points values where to evaluate the Morison element contribution (local coordinates referred to the COG)
    cut::CusTensor<cusfloat>                    _field_points_g     ;           // Field points values where to evaluate the Morison element contribution (global coordinates)
    std::size_t                                 _field_points_np    = 0;        // Number of field points
    Input*                                      _input              = nullptr;  // Pointer to input data
    cusfloat                                    _kc_length          = 0.0;      // Characteristic length used for the Keulegan-Carpenter number calculation
    RadDiffData<cuscomplex, RDDMorisonConfig>*  _rdd_morison        ;           // Radiation and diffraction data for Morison element calculations
    MorisonCoeffCurve*                          _cd_axial           = nullptr;  // Drag coefficient for the axial component of the Morison element
    MorisonCoeffCurve*                          _cd_trans           = nullptr;  // Drag coefficient for the transverse component of the Morison element
    MorisonCoeffCurve*                          _cd_vert            = nullptr;  // Drag coefficient for the vertical component of the Morison element
    cut::CusTensor<cusfloat>                    _x_axis_l           ;           // Local x axis of the Morison element
    cut::CusTensor<cusfloat>                    _y_axis_l           ;           // Local y axis of the Morison element
    cut::CusTensor<cusfloat>                    _z_axis_l           ;           // Local z axis of the Morison element


    /*** Define private class methods ***/
    void    _calculate_field_points( 
                                        MorisonElementDef*  morison_def_
                                    );

    void    _initialize_element( 
                                    MpiConfig*          mpi_config_,
                                    Input*              input_,
                                    MorisonElementDef*  morison_def_
                                );


public:
    /*** Class constructors declaration ***/
    MorisonElement( ) = default;

    MorisonElement( 
                        MpiConfig*          mpi_config_,
                        Input*              input_,
                        MorisonElementDef*  morison_def_
                    );

    MorisonElement(const MorisonElement&) = delete;
    MorisonElement& operator=(const MorisonElement&) = delete;

    MorisonElement(MorisonElement&&) = default;
    MorisonElement& operator=(MorisonElement&&) = default;

    ~MorisonElement( );

    /*** Class public methods declaration ***/
    void calculate_hydrodynamic_forces( 
                                            cusfloat                    ang_freq,
                                            cuscomplex*                 raos,
                                            cut::CusTensor<cuscomplex>& inertial_force,
                                            cut::CusTensor<cuscomplex>& drag_force
                                        ) const;

};