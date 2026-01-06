
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


/***********************************************************/
/**************** Module Auxiliar Macros *******************/
/***********************************************************/
#define _DELETE_TENSOR_FIELD( tensor_ptr )              \
    if ( tensor_ptr != nullptr )                        \
    {                                                   \
        delete tensor_ptr;                              \
        tensor_ptr = nullptr;                           \
    }                                                   \


/***********************************************************/
/**************** Define RadDiffData class *****************/
/***********************************************************/
template<int mode_f, int mode_dfdn, int mode_dfdc>
RadDiffData<mode_f, mode_dfdn, mode_dfdc>::RadDiffData( 
                                                            std::size_t     field_points_np_,
                                                            std::size_t     dofs_np,
                                                            std::size_t     headings_np
                                                        )
{
    // Store number of field points
    this->field_points_np   = field_points_np_;

    // Allocate memory on heap for field points
    this->field_points      = new cut::CusTensor<cusfloat>( { field_points_np_, 3 } );

    // Allocate memory on heap for potential fields
    STATIC_COND( ONLY_FCN,   this->pot_incident    = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCN,   this->pot_raddiff     = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCN,   this->pot_total       = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )

    // Allocate memory on heap for normal velocity derivative fields
    STATIC_COND( ONLY_FCNDN, this->vel_dn_incident = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCNDN, this->vel_dn_raddiff  = new cut::CusTensor<cuscomplex>( { field_points_np_, dofs_np * headings_np } ); )
    STATIC_COND( ONLY_FCNDN, this->vel_dn_total    = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )

    // Allocate memory on heap for velocity components fields
    STATIC_COND( ONLY_FCNDC, this->vel_x_incident  = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCNDC, this->vel_y_incident  = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCNDC, this->vel_z_incident  = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCNDC, this->vel_x_raddiff   = new cut::CusTensor<cuscomplex>( { field_points_np_, dofs_np * headings_np } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_y_raddiff   = new cut::CusTensor<cuscomplex>( { field_points_np_, dofs_np * headings_np } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_z_raddiff   = new cut::CusTensor<cuscomplex>( { field_points_np_, dofs_np * headings_np } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_x_total     = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCNDC, this->vel_y_total     = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )
    STATIC_COND( ONLY_FCNDC, this->vel_z_total     = new cut::CusTensor<cuscomplex>( { field_points_np_, headings_np } );           )

    // Set flag to indicate that memory is allocated on heap
    this->_is_heap        = true;

}


template<int mode_f, int mode_dfdn, int mode_dfdc>
RadDiffData<mode_f, mode_dfdn, mode_dfdc>::~RadDiffData( )
{
    // Deallocate memory only if it was allocated on heap
    if ( this->_is_heap )
    {
        // Delete field points memory
        _DELETE_TENSOR_FIELD( this->field_points )

        // Delete potential fields memory
        STATIC_COND( ONLY_FCN,   _DELETE_TENSOR_FIELD( this->pot_incident )    )
        STATIC_COND( ONLY_FCN,   _DELETE_TENSOR_FIELD( this->pot_raddiff  )    )
        STATIC_COND( ONLY_FCN,   _DELETE_TENSOR_FIELD( this->pot_total    )    )

        // Delete normal velocity derivative fields memory
        STATIC_COND( ONLY_FCNDN, _DELETE_TENSOR_FIELD( this->vel_dn_incident ) )
        STATIC_COND( ONLY_FCNDN, _DELETE_TENSOR_FIELD( this->vel_dn_raddiff  ) )
        STATIC_COND( ONLY_FCNDN, _DELETE_TENSOR_FIELD( this->vel_dn_total    ) )

        // Delete velocity fields components memory
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_x_incident )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_y_incident )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_z_incident )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_x_raddiff  )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_y_raddiff  )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_z_raddiff  )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_x_total    )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_y_total    )  )
        STATIC_COND( ONLY_FCNDC, _DELETE_TENSOR_FIELD( this->vel_z_total    )  )

    }
}