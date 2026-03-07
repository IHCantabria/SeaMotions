
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

// Include general usage libraries
#include <filesystem>
#include <string>

// Include local modules
#include "../math/data_fit/interp1d.hpp"
#include "../math/custensor/custensor.hpp"


struct MorisonCoeffCurve
{
private:
    /*** Define private class attributes ***/
    Interp1D<cusfloat>*         _interp   = nullptr;  // Added mass coefficient
    cut::CusTensor<cusfloat>    _x_data;                // X data of the curve
    cut::CusTensor<cusfloat>    _y_data;                // Y data of the curve


    /*** Define private class methods ***/
    void _parse_input_file( 
                                std::string foname,
                                std::string finame
                            );

public:
    /*** Define public class attributes ***/

    /*** Define public class methods ***/
    MorisonCoeffCurve( ) = default;

    MorisonCoeffCurve( 
                        std::string foname,
                        std::string finame
                    );

    MorisonCoeffCurve( const MorisonCoeffCurve& ) = delete;
    MorisonCoeffCurve& operator=( const MorisonCoeffCurve& ) = delete;

    MorisonCoeffCurve( MorisonCoeffCurve&& ) = default;
    MorisonCoeffCurve& operator=( MorisonCoeffCurve&& ) = default;

    ~MorisonCoeffCurve( );

    /*** Declare public methods ***/
    cusfloat eval( cusfloat x ) const;

};