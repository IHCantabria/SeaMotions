
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

// Include general usage libraries
#include <vector>

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
template<int mode_comp, int mode_f, int mode_dfdn, int mode_dfdc>
void    PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>::_allocate_memory( 
                                                                                std::size_t field_points_np_,
                                                                                std::size_t freqs_np_,
                                                                                std::size_t headings_np_,
                                                                                std::size_t dofs_np_
                                                                            )
{
    // Set number of field points
    this->dofs_np           = dofs_np_;
    this->field_points_np   = field_points_np_;
    this->freqs_np          = freqs_np_;
    this->headings_np       = headings_np_;

    // Define vector shape sizes for fields
    std::vector<std::size_t> shape_fp( { freqs_np_, headings_np_, field_points_np_ } );

    // Allocate memory on heap for field points
    this->field_points      = new cut::CusTensor<cusfloat>( { field_points_np_, 3 } );

    // Allocate memory on heap for potential fields
    STATIC_COND( ONLY_FCN && USE_COMP,   this->pot_incident     = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCN && USE_COMP,   this->pot_rad          = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCN && USE_COMP,   this->pot_diff         = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCN            ,   this->pot_total        = new cut::CusTensor<cuscomplex>( shape_fp ); )

    // Allocate memory on heap for pressure field
    STATIC_COND( ONLY_FCN,               this->press_total      = new cut::CusTensor<cuscomplex>( shape_fp ); )

    // Allocate memory on heap for normal velocity derivative fields
    STATIC_COND( ONLY_FCNDN && USE_COMP, this->vel_dn_incident  = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDN && USE_COMP, this->vel_dn_rad       = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDN && USE_COMP, this->vel_dn_diff      = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDN,             this->vel_dn_total     = new cut::CusTensor<cuscomplex>( shape_fp ); )

    // Allocate memory on heap for velocity components fields
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_x_incident   = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_y_incident   = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_z_incident   = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_x_rad        = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_y_rad        = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_z_rad        = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_x_diff       = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_y_diff       = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP, this->vel_z_diff       = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC,             this->vel_x_total      = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC,             this->vel_y_total      = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCNDC,             this->vel_z_total      = new cut::CusTensor<cuscomplex>( shape_fp ); )

    // Allocate memory on heap for wave elevation fields
    STATIC_COND( ONLY_FCN,               this->wev_total        = new cut::CusTensor<cuscomplex>( shape_fp ); )
    STATIC_COND( ONLY_FCN,               this->wev_rel_total    = new cut::CusTensor<cuscomplex>( shape_fp ); )

    // Set flag to indicate that memory is allocated on heap
    this->_is_heap        = true;

}


template<int mode_comp, int mode_f, int mode_dfdn, int mode_dfdc>
void    PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>::_load_field_points( 
                                                                                PanelGeom*  panel_geom_,
                                                                                bool        use_waterline_
                                                                            )
{
    if ( use_waterline_ )
    {
        int         count_lines = 0;
        cusfloat    tv[3];      clear_vector( 3, tv );
        GaussPoints gp( panel_geom_->gauss_points_np );

        // Create direction vector
        tv[0]   = ( panel_geom_->x_wl[1] - panel_geom_->x_wl[0] ) / panel_geom_->len_wl;
        tv[1]   = ( panel_geom_->y_wl[1] - panel_geom_->y_wl[0] ) / panel_geom_->len_wl;

        for ( std::size_t i=0; i<panel_geom_->gauss_points_np; i++ )
        {
            this->field_points[3*i]     = tv[0] * ( gp.roots[i] + 1.0 ) / 2.0 * panel_geom_->len_wl + panel_geom_->x_wl[0];
            this->field_points[3*i+1]   = tv[1] * ( gp.roots[i] + 1.0 ) / 2.0 * panel_geom_->len_wl + panel_geom_->y_wl[0];
            this->field_points[3*i+2]   = 0.0;
        }

    }
    else
    {
        // Loop over gauss points defined in PanelGeom to 
        // load field points coordinates
        for ( std::size_t i = 0; i < pow2s( panel_geom_->gauss_points_np ); ++i )
        {
            this->field_points[3*i + 0] = panel_geom_->gauss_points_global_x[i];
            this->field_points[3*i + 1] = panel_geom_->gauss_points_global_y[i];
            this->field_points[3*i + 2] = panel_geom_->gauss_points_global_z[i];
        }
    }
}

                                        
template<int mode_comp, int mode_f, int mode_dfdn, int mode_dfdc>
void PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>::clear_data( void )
{
    // Create auxiliary variables to have a 
    // clear implementation
    std::size_t hfnp = this->freqs_np * this->headings_np * this->field_points_np;

    // Clear potential fields data
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->pot_incident[i]      = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->pot_rad[i]           = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->pot_diff[i]          = 0.0; ) )
    STATIC_COND( ONLY_FCN,                  LOOP_DEF( hfnp,  this->pot_total[i]         = 0.0; ) )

    // Clear normal velocity derivative fields data
    STATIC_COND( ONLY_FCNDN && USE_COMP,    LOOP_DEF( hfnp,  this->vel_dn_incident[i]   = 0.0; ) )
    STATIC_COND( ONLY_FCNDN && USE_COMP,    LOOP_DEF( hfnp,  this->vel_dn_raddiff[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDN && USE_COMP,    LOOP_DEF( hfnp,  this->vel_dn_raddiff[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDN,                LOOP_DEF( hfnp,  this->vel_dn_total[i]      = 0.0; ) )

    // Clear velocity components fields data
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_x_incident[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_y_incident[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_z_incident[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_x_rad[i]         = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_y_rad[i]         = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_z_rad[i]         = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_x_diff[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_y_diff[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_COMP,    LOOP_DEF( hfnp,  this->vel_z_diff[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCNDC,                LOOP_DEF( hfnp,  this->vel_x_total[i]       = 0.0; ) )
    STATIC_COND( ONLY_FCNDC,                LOOP_DEF( hfnp,  this->vel_y_total[i]       = 0.0; ) )
    STATIC_COND( ONLY_FCNDC,                LOOP_DEF( hfnp,  this->vel_z_total[i]       = 0.0; ) )
}

template<int mode_comp, int mode_f, int mode_dfdn, int mode_dfdc>
PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>::PanelData( 
                                                        std::size_t field_points_np_,
                                                        std::size_t freqs_np_,
                                                        std::size_t headings_np_,
                                                        std::size_t dofs_np_
                                                    )
{
    // Allocate space for fields data
    this->_allocate_memory( 
                                field_points_np_,
                                freqs_np_,
                                headings_np_,
                                dofs_np_ 
                            );
}


template<int mode_comp, int mode_f, int mode_dfdn, int mode_dfdc>
PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>::PanelData( 
                                                                PanelGeom*  panel_geom_,
                                                                std::size_t body_id_,
                                                                std::size_t freqs_np_,
                                                                std::size_t headings_np_,
                                                                std::size_t dofs_np_,
                                                                bool        use_waterline_
                                                            )
{
    // Storage input arguments
    this->panel_geom        = panel_geom_;
    this->body_id           = body_id_;

    // Allocate space for fields data
    std::size_t gp_np       = panel_geom_->gauss_points_np;
    std::size_t fp_np       = ( use_waterline_) ? gp_np : pow2s(  gp_np );

    this->_allocate_memory( 
                                fp_np,
                                freqs_np_,
                                headings_np_,
                                dofs_np_ 
                            );

    // Load field points from PanelGeom
    this->_load_field_points( 
                                this->_panel_geom,
                                use_waterline_
                            );

}


template<int mode_comp, int mode_f, int mode_dfdn, int mode_dfdc>
PanelData<mode_comp, mode_f, mode_dfdn, mode_dfdc>::~PanelData( )
{
    // Deallocate memory only if it was allocated on heap
    if ( this->_is_heap )
    {
        // Delete field points memory
        _DELETE_TENSOR_FIELD( this->field_points )

        // Delete potential fields memory
        STATIC_COND( ONLY_FCN && USE_COMP,      _DELETE_TENSOR_FIELD( this->pot_incident    ) )
        STATIC_COND( ONLY_FCN && USE_COMP,      _DELETE_TENSOR_FIELD( this->pot_rad         ) )
        STATIC_COND( ONLY_FCN && USE_COMP,      _DELETE_TENSOR_FIELD( this->pot_diff        ) )
        STATIC_COND( ONLY_FCN,                  _DELETE_TENSOR_FIELD( this->pot_total       ) )

        // Delete normal velocity derivative fields memory
        STATIC_COND( ONLY_FCNDN && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_dn_incident ) )
        STATIC_COND( ONLY_FCNDN && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_dn_rad      ) )
        STATIC_COND( ONLY_FCNDN && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_dn_diff     ) )
        STATIC_COND( ONLY_FCNDN,                _DELETE_TENSOR_FIELD( this->vel_dn_total    ) )

        // Delete velocity fields components memory
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_x_incident  ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_y_incident  ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_z_incident  ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_x_rad       ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_y_rad       ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_z_rad       ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_x_diff      ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_y_diff      ) )
        STATIC_COND( ONLY_FCNDC && USE_COMP,    _DELETE_TENSOR_FIELD( this->vel_z_diff      ) )
        STATIC_COND( ONLY_FCNDC,                _DELETE_TENSOR_FIELD( this->vel_x_total     ) )
        STATIC_COND( ONLY_FCNDC,                _DELETE_TENSOR_FIELD( this->vel_y_total     ) )
        STATIC_COND( ONLY_FCNDC,                _DELETE_TENSOR_FIELD( this->vel_z_total     ) )

    }
}