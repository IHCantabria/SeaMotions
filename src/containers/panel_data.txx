
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
#include "panel_data.hpp"
#include "../static_tools.hpp"


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
void PanelData<mode_f, mode_dfdn, mode_dfdc>::clear_data( void )
{
    // Create auxiliary variables to have a 
    // clear implementation
    std::size_t dhfnp    = this->dofs_np * this->headings_np * this->field_points_np;
    std::size_t hfnp     = this->headings_np * this->field_points_np;

    // Clear potential fields data
    STATIC_COND( ONLY_FCN,   LOOP_DEF( hfnp,  this->pot_incident[i]     = 0.0; ) )
    STATIC_COND( ONLY_FCN,   LOOP_DEF( dhfnp, this->pot_raddiff[i]      = 0.0; ) )
    STATIC_COND( ONLY_FCN,   LOOP_DEF( hfnp,  this->pot_total[i]        = 0.0; ) )

    // Clear normal velocity derivative fields data
    STATIC_COND( ONLY_FCNDN, LOOP_DEF( hfnp,  this->vel_dn_incident[i]  = 0.0; ) )
    STATIC_COND( ONLY_FCNDN, LOOP_DEF( dhfnp, this->vel_dn_raddiff[i]   = 0.0; ) )
    STATIC_COND( ONLY_FCNDN, LOOP_DEF( hfnp,  this->vel_dn_total[i]     = 0.0; ) )

    // Clear velocity components fields data
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( hfnp,  this->vel_x_incident[i]   = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( hfnp,  this->vel_y_incident[i]   = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( hfnp,  this->vel_z_incident[i]   = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( dhfnp, this->vel_x_raddiff[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( dhfnp, this->vel_y_raddiff[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( dhfnp, this->vel_z_raddiff[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( hfnp,  this->vel_x_total[i]      = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( hfnp,  this->vel_y_total[i]      = 0.0; ) )
    STATIC_COND( ONLY_FCNDC, LOOP_DEF( hfnp,  this->vel_z_total[i]      = 0.0; ) )
}

template<int mode_f, int mode_dfdn, int mode_dfdc>
PanelData<mode_f, mode_dfdn, mode_dfdc>::PanelData( 
                                                        std::size_t field_points_np_,
                                                        std::size_t headings_np_,
                                                        std::size_t dofs_np_
                                                    )
{
    // Set number of field points
    this->dofs_np           = dofs_np_;
    this->field_points_np   = field_points_np_;
    this->headings_np       = headings_np_;

    // Allocate memory on heap for field points
    this->field_points      = new cut::CusTensor<cusfloat>( { field_points_np_, 3 } );

    // Allocate memory on heap for potential fields
    STATIC_COND( ONLY_FCN,   this->pot_incident    = new cut::CusTensor<cuscomplex>( { headings_np ,          field_points_np_ } ); )
    STATIC_COND( ONLY_FCN,   this->pot_raddiff     = new cut::CusTensor<cuscomplex>( { dofs_np + headings_np, field_points_np_ } ); )
    STATIC_COND( ONLY_FCN,   this->pot_total       = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )

    // Allocate memory on heap for normal velocity derivative fields
    STATIC_COND( ONLY_FCNDN, this->vel_dn_incident = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDN, this->vel_dn_raddiff  = new cut::CusTensor<cuscomplex>( { dofs_np + headings_np, field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDN, this->vel_dn_total    = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )

    // Allocate memory on heap for velocity components fields
    STATIC_COND( ONLY_FCNDC, this->vel_x_incident  = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_y_incident  = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_z_incident  = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_x_raddiff   = new cut::CusTensor<cuscomplex>( { dofs_np + headings_np, field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_y_raddiff   = new cut::CusTensor<cuscomplex>( { dofs_np + headings_np, field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_z_raddiff   = new cut::CusTensor<cuscomplex>( { dofs_np + headings_np, field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_x_total     = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_y_total     = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )
    STATIC_COND( ONLY_FCNDC, this->vel_z_total     = new cut::CusTensor<cuscomplex>( { headings_np,           field_points_np_ } ); )

    // Set flag to indicate that memory is allocated on heap
    this->_is_heap        = true;
}


template<int mode_f, int mode_dfdn, int mode_dfdc>
PanelData<mode_f, mode_dfdn, mode_dfdc>::~PanelData( )
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