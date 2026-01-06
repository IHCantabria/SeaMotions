
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

// Include local modules
#include "rad_diff_data.hpp"


RadDiffData::RadDiffData( 
                                std::size_t     field_points_np_,
                                std::size_t     dofs_np,
                                std::size_t     headings_np
                            )
{
    // Store number of field points
    this->field_points_np = field_points_np_;

    // Allocate memory on heap for field points
    this->field_points    = new cut::CusTensor<cusfloat>( { field_points_np_, 3 } );

    // Allocate memory on heap for incident field
    this->incident        = new cut::CusTensor<cuscomplex>( {field_points_np_, headings_np } );

    // Allocate memory on heap for radiated + diffracted field
    this->raddiff         = new cut::CusTensor<cuscomplex>( { field_points_np_, dofs_np * headings_np } );

    // Allocate memory on heap for total field
    this->total           = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );

    // Set flag to indicate that memory is allocated on heap
    this->_is_heap        = true;

}


RadDiffData::~RadDiffData( )
{
    // Deallocate memory only if it was allocated on heap
    if ( this->_is_heap )
    {
        // Delete field points memory
        if ( this->field_points != nullptr )
        {
            delete this->field_points;
            this->field_points = nullptr;
        }

        // Delete incident field memory
        if ( this->incident != nullptr )
        {
            delete this->incident;
            this->incident = nullptr;
        }

        // Delete radiated + diffracted field memory
        if ( this->raddiff != nullptr )
        {
            delete this->raddiff;
            this->raddiff = nullptr;
        }

        // Delete total field memory
        if ( this->total != nullptr )
        {
            delete this->total;
            this->total = nullptr;
        }
    }
}