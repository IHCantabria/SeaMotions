
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
#include "../math/gauss.hpp"
#include "panel_data.hpp"
#include "../static_tools.hpp"


/***********************************************************/
/**************** Define RadDiffData class *****************/
/***********************************************************/
template<typename T, typename Config>
void    PanelData<T, Config>::_allocate_memory( 
                                                std::size_t field_points_np_,
                                                std::size_t freqs_np_,
                                                std::size_t headings_np_,
                                                std::size_t dofs_np_
                                            )
{
    // Set number of field points
    this->dofs_np           = dofs_np_;
    this->field_points_np   = field_points_np_;
    this->freqs_np          = ( store_freqs == ON ) ? freqs_np_ : 1;
    this->headings_np       = headings_np_;

    // Define vector shape sizes for fields
    std::size_t field_len   = this->freqs_np * this->headings_np * this->field_points_np;
    std::size_t mode_len    = this->freqs_np * ( dofs_np_ + headings_np_ ) * this->field_points_np;

    // Allocate memory on heap for field points
    this->field_points.resize( 3 * field_points_np_ );

    // Allocate memory on heap for potential fields
    STATIC_COND( ONLY_FCN && USE_COMP,   this->pot_incident.resize( field_len );    )
    STATIC_COND( ONLY_FCN && USE_COMP,   this->pot_rad.resize( field_len );         )
    STATIC_COND( ONLY_FCN && USE_COMP,   this->pot_diff.resize( field_len );        )
    STATIC_COND( ONLY_FCN,               this->pot_total.resize( field_len );       )
    STATIC_COND( ONLY_FCN,               this->pot_total_2.resize( field_len );     )
    STATIC_COND( ONLY_FCN && USE_MODES,  this->pot_modes.resize( mode_len );        )

    // Allocate memory on heap for pressure field
    STATIC_COND( ONLY_FCN && USE_COMP,   this->press_incident.resize( field_len );  )
    STATIC_COND( ONLY_FCN && USE_COMP,   this->press_rad.resize( field_len );       )
    STATIC_COND( ONLY_FCN && USE_COMP,   this->press_diff.resize( field_len );      )
    STATIC_COND( ONLY_FCN,               this->press_total.resize( field_len );     )

    // Allocate memory on heap for normal velocity derivative fields
    STATIC_COND( ONLY_FCNDN && USE_COMP, this->vel_dn_incident.resize( field_len ); )
    STATIC_COND( ONLY_FCNDN && USE_COMP, this->vel_dn_rad.resize( field_len );      )
    STATIC_COND( ONLY_FCNDN && USE_COMP, this->vel_dn_diff.resize( field_len );     )
    STATIC_COND( ONLY_FCNDN,             this->vel_dn_total.resize( field_len );    )

    // Allocate memory on heap for velocity components fields
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_x_incident.resize( field_len ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_y_incident.resize( field_len ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_z_incident.resize( field_len ); )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_x_rad.resize( field_len );      )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_y_rad.resize( field_len );      )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_z_rad.resize( field_len );      )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_x_diff.resize( field_len );     )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_y_diff.resize( field_len );     )
    STATIC_COND( ONLY_FCNDC && USE_COMP,  this->vel_z_diff.resize( field_len );     )
    STATIC_COND( ONLY_FCNDC,              this->vel_x_total.resize( field_len );    )
    STATIC_COND( ONLY_FCNDC,              this->vel_y_total.resize( field_len );    )
    STATIC_COND( ONLY_FCNDC,              this->vel_z_total.resize( field_len );    )
    STATIC_COND( ONLY_FCNDC && USE_MODES, this->vel_x_modes.resize( mode_len );     )
    STATIC_COND( ONLY_FCNDC && USE_MODES, this->vel_y_modes.resize( mode_len );     )
    STATIC_COND( ONLY_FCNDC && USE_MODES, this->vel_z_modes.resize( mode_len );     )

    // Allocate memory on heap for wave elevation fields
    STATIC_COND( ONLY_FCN,               this->wev_total.resize( field_len );       )
    STATIC_COND( ONLY_FCN,               this->wev_rel_total.resize( field_len );   )

    // Set flag to indicate that memory is allocated on heap
    this->_is_heap        = true;

}


template<typename T, typename Config>
template<FieldTypeE field_type, FieldComponentE field_component>
const cut::CusTensor<T>* PanelData<T, Config>::get_field_data( void ) const
{
    if constexpr ( field_type == FieldTypeE::POTENTIAL )
    {
        if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::INCIDENT ) )
            return &(this->pot_incident);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::RADIATED ) )
            return &(this->pot_rad);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::SCATTERED ) )
            return &(this->pot_diff);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->pot_total);
    }
    else if constexpr ( field_type == FieldTypeE::PRESSURE )
    {
        if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::INCIDENT ) )
            return &(this->press_incident);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::RADIATED ) )
            return &(this->press_rad);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::SCATTERED ) )
            return &(this->press_diff);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->press_total);
    }
    else if constexpr ( field_type == FieldTypeE::VELOCITY_DN )
    {
        if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::INCIDENT ) )
            return &(this->vel_dn_incident);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::RADIATED ) )
            return &(this->vel_dn_rad);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::SCATTERED ) )
            return &(this->vel_dn_diff);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->vel_dn_total);
    }
    else if constexpr ( field_type == FieldTypeE::VELOCITY_X )
    {
        if constexpr ( ( mode_comp == ON ) && (  field_component == FieldComponentE::INCIDENT ) )
            return &(this->vel_x_incident);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::RADIATED ) )
            return &(this->vel_x_rad);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::SCATTERED ) )
            return &(this->vel_x_diff);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->vel_x_total);
    }
    else if constexpr ( field_type == FieldTypeE::VELOCITY_Y )
    {
        if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::INCIDENT ) )
            return &(this->vel_y_incident);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::RADIATED ) )
            return &(this->vel_y_rad);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::SCATTERED ) )
            return &(this->vel_y_diff);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->vel_y_total);
    }
    else if constexpr ( field_type == FieldTypeE::VELOCITY_Z )
    {
        if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::INCIDENT ) )
            return &(this->vel_z_incident);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::RADIATED ) )
            return &(this->vel_z_rad);
        else if constexpr ( ( mode_comp == ON ) && ( field_component == FieldComponentE::SCATTERED ) )
            return &(this->vel_z_diff);
        else if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->vel_z_total);
    }
    else if constexpr ( field_type == FieldTypeE::WAVE_ELEVATION )
    {
        if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->wev_total);
    }
    else if constexpr ( field_type == FieldTypeE::RELATIVE_WAVE_ELEVATION )
    {
        if constexpr ( ( mode_f == ON ) && ( field_component == FieldComponentE::TOTAL ) )
            return &(this->wev_rel_total);
    }
    return nullptr;
}

template<typename T, typename Config>
void    PanelData<T, Config>::_load_field_points( 
                                                    PanelGeom*  panel_geom_,
                                                    bool        use_waterline_
                                                )
{
    if ( use_waterline_ )
    {
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

                                        
template<typename T, typename Config>
void PanelData<T, Config>::clear_data( void )
{
    // Create auxiliary variables to have a 
    // clear implementation
    std::size_t hfnp = this->freqs_np * this->headings_np * this->field_points_np;
    std::size_t mfnp = this->freqs_np * ( this->dofs_np + this->headings_np ) * this->field_points_np;

    // Clear potential fields data
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->pot_incident[i]      = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->pot_rad[i]           = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->pot_diff[i]          = 0.0; ) )
    STATIC_COND( ONLY_FCN,                  LOOP_DEF( hfnp,  this->pot_total[i]         = 0.0; ) )
    STATIC_COND( ONLY_FCN,                  LOOP_DEF( hfnp,  this->pot_total_2[i]       = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_MODES,     LOOP_DEF( mfnp,  this->pot_modes[i]         = 0.0; ) )

    // Clear pressure fields data
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->press_incident[i]    = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->press_rad[i]         = 0.0; ) )
    STATIC_COND( ONLY_FCN && USE_COMP,      LOOP_DEF( hfnp,  this->press_diff[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCN,                  LOOP_DEF( hfnp,  this->press_total[i]       = 0.0; ) )

    // Clear normal velocity derivative fields data
    STATIC_COND( ONLY_FCNDN && USE_COMP,    LOOP_DEF( hfnp,  this->vel_dn_incident[i]   = 0.0; ) )
    STATIC_COND( ONLY_FCNDN && USE_COMP,    LOOP_DEF( hfnp,  this->vel_dn_rad[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCNDN && USE_COMP,    LOOP_DEF( hfnp,  this->vel_dn_diff[i]       = 0.0; ) )
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
    STATIC_COND( ONLY_FCNDC && USE_MODES,   LOOP_DEF( mfnp, this->vel_x_modes[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_MODES,   LOOP_DEF( mfnp, this->vel_y_modes[i]        = 0.0; ) )
    STATIC_COND( ONLY_FCNDC && USE_MODES,   LOOP_DEF( mfnp, this->vel_z_modes[i]        = 0.0; ) )

    // Clear wave elevation fields data
    STATIC_COND( ONLY_FCN,                  LOOP_DEF( hfnp,  this->wev_total[i]         = 0.0; ) )
    STATIC_COND( ONLY_FCN,                  LOOP_DEF( hfnp,  this->wev_rel_total[i]     = 0.0; ) )

}

template<typename T, typename Config>
PanelData<T, Config>::PanelData( 
                                    std::size_t field_points_np_,
                                    cusfloat*   field_points_,
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

    // Copy field points coordinates from input array
    for (  std::size_t i = 0; i < 3 * field_points_np_; ++i )
    {
        this->field_points[i] = field_points_[i];
    }
}


template<typename T, typename Config>
PanelData<T, Config>::PanelData( 
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
                                this->panel_geom,
                                use_waterline_
                            );

}