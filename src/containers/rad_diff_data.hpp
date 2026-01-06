
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
#include "../math/custensor/custensor.hpp"


struct RadDiffData
{
private:
    // Declare private variables
    bool                            _is_heap        = false; // Flag to indicate if memory is allocated on heap

public:
    // Declare public variables
    cut::CusTensor<cusfloat>*       field_points    = nullptr; // Store field points coordinates
    std::size_t                     field_points_np = 0;       // Number of field points
    cut::CusTensor<cuscomplex>*     incident        = nullptr; // Store wave incident field value [field_points, headings_np]
    cut::CusTensor<cuscomplex>*     raddiff         = nullptr; // Store radiated diffracted field values [field_points, dofs_np · headings_np]
    cut::CusTensor<cuscomplex>*     total           = nullptr; // Store total field value [field_points, headings_np]. Total field is composed of: incident + sum( epsi · radiation ) + diffraction.

    /* Declare class constructors and destructor */
    RadDiffData( ) = default;

    RadDiffData( 
                    std::size_t     field_points_np,
                    std::size_t     dofs_np,
                    std::size_t     headings_np
                );

    ~RadDiffData( );
    
};
