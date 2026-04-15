
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
#include "qtf.hpp"


inline void     _calculate_so_rotation_mat( 
                                            cuscomplex* rots_i,
                                            cuscomplex* rots_j,
                                            QTFTypeE    qtf_type,
                                            cuscomplex* H
                                        )
{
    // Get local variables
    cuscomplex  ai[3]   = { rots_i[0], rots_i[1], rots_i[2] };
    cuscomplex  aj[3]   = { rots_j[0], rots_j[1], rots_j[2] };

    // Conjugate j values if difference frequency case is required
    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        aj[0] = std::conj( aj[0] );
        aj[1] = std::conj( aj[1] );
        aj[2] = std::conj( aj[2] );
    }

    // Define rotation matrix
    H[0] = - 0.5 * ( ai[1] * aj[1] + ai[2] * aj[2] );
    H[1] = 0.0;
    H[2] = 0.0;
    H[3] = + 0.5 * ( ai[0] * aj[1] + ai[1] * aj[0] );
    H[4] = - 0.5 * ( ai[0] * aj[0] + ai[2] * aj[2] );
    H[5] = 0.0;
    H[6] = + 0.5 * ( ai[0] * aj[2] + ai[2] * aj[0] );
    H[7] = + 0.5 * ( ai[1] * aj[2] + ai[2] * aj[1] );
    H[8] = - 0.5 * ( ai[0] * aj[0] + ai[1] * aj[1] );
}


template<typename Config>
inline  void    so_potential_qb_rhs(
                                            QTFTypeE                        qtf_type,
                                            Input*                          input,
                                            std::size_t                     freq_i,
                                            std::size_t                     freq_j,
                                            cuscomplex*                     raos_hist,
                                            std::size_t                     body_id,
                                            cusfloat*                       cog,
                                            PanelData<cuscomplex, Config>*  panel_data,
                                            cuscomplex*                     G,
                                            cuscomplex*                     dG_dx,
                                            cuscomplex*                     dG_dy,
                                            cuscomplex*                     dG_dz,
                                            WaveDispersionSO*               wd,
                                            bool                            mode_wl,
                                            cuscomplex*                     qb_rhs
                                )
{
    // Get local variables
    std::size_t bodies_np       = static_cast<input->bodies_np>;
    std::size_t dofs_np         = static_cast<input->dofs_np>;
    std::size_t heads_np        = static_cast<input->heads_np>;
    cuscomplex  dvdn            = 0.0;
    cusfloat    field_point[3]  = { 0.0, 0.0, 0.0 };
    cuscomplex  fp_cog[3]       = { 0.0, 0.0, 0.0 };
    cuscomplex  fp_so_i[3]      = { 0.0, 0.0, 0.0 };
    cuscomplex  fp_so_j[3]      = { 0.0, 0.0, 0.0 };
    std::size_t index           = 0;
    cusfloat*   nvec            = panel_data->panel_geom->normal_vec;
    std::size_t pd_idx          = 0;
    std::size_t pd_jdx          = 0;
    cuscomplex* rao_i           = &raos_freq_i[ freq_i * heads_np * bodies_np * dofs_np ];
    cuscomplex* rao_j           = &raos_freq_j[ freq_j * heads_np * bodies_np * dofs_np ];
    std::size_t rao_i_idx       = 0;
    std::size_t rao_j_idx       = 0;
    cuscomplex  so_rotmat[9]    ;
    cusfloat    tv[3]           = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_i[3]        = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_j[3]        = { 0.0, 0.0, 0.0 };
    cuscomplex  v0[3]           = { 0.0, 0.0, 0.0 };
    cuscomplex  v1[3]           = { 0.0, 0.0, 0.0 };
    cuscomplex  wave_dx         = 0.0;
    cuscomplex  wave_dy         = 0.0;
    cuscomplex  wave_dz         = 0.0;

    // Define QTF sign based on the type of QTF term
    int qtf_sign                = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? -1 : 1;

    // Get angular frequencies
    cusfloat wi                 = input->ang_freqs[freq_i];
    cusfloat wj                 = input->ang_freqs[freq_j];
    cusfloat wi_wj              = wi + qtf_sign * wj;

    // Clear rhs vector to sum up all the constributions on it
    std::size_t size_rhs        = pow2s( heads_np ) * panel_data->field_points_np;
    clear_vector( size_rhs, qb_rhs );

    // Loop over headings to get all the possible combinations
    for ( std::size_t ih=0; ih<heads_np; ih++ )
    {
        for ( std::size_t jh=0; jh<heads_np; jh++ )
        {
            // Get RAO values for the current frequency and heading combination
            rao_i_idx   = freq_i * heads_np * bodies_np * dofs_np + ih * bodies_np * dofs_np + body_id * dofs_np;
            rao_j_idx   = freq_j * heads_np * bodies_np * dofs_np + jh * bodies_np * dofs_np + body_id * dofs_np;
            rao_i       = &raos_freq_i[ rao_i_idx ];
            rao_j       = &raos_freq_j[ rao_j_idx ];

            // Update wave dispersion relation for second order wave propagation
            wd->set_new_data( 
                                wi, 
                                wj, 
                                panel_data->head[ih], 
                                panel_data->head[jh] 
                            );

            for ( std::size_t fpi=0; fpi<panel_data->field_points_np; fpi++ )
            {
                index           = ( ih * heads_np + jh ) * panel_data->field_points_np + fpi;
                // Get current field point coordinates
                field_point[0]  = panel_data->field_points[3*fpi+0];
                field_point[1]  = panel_data->field_points[3*fpi+1];
                field_point[2]  = panel_data->field_points[3*fpi+2];

                // Get current field point velocities for both frequencies combinations
                pd_idx          = freq_i * heads_np * panel_data->field_points_np;
                vel_i[0]        = panel_data->vel_x_total[pd_idx + fpi];
                vel_i[1]        = panel_data->vel_y_total[pd_idx + fpi];
                vel_i[2]        = panel_data->vel_z_total[pd_idx + fpi];
                
                pd_jdx          = freq_j * heads_np * panel_data->field_points_np;
                vel_j[0]        = panel_data->vel_x_total[pd_jdx + fpi];
                vel_j[1]        = panel_data->vel_y_total[pd_jdx + fpi];
                vel_j[2]        = panel_data->vel_z_total[pd_jdx + fpi];

                // Define field point coordinates with respect to body cog
                fp_cog[0]       = field_point[0] - cog[0];
                fp_cog[1]       = field_point[1] - cog[1];
                fp_cog[2]       = field_point[2] - cog[2];

                clear_vector( 3, fp_so_i );
                cross( rao_i, fp_cog, fp_so_i );
                sv_add( 3, fp_so_i, rao_i, fp_so_i );

                clear_vector( 3, fp_so_j );
                cross( rao_j, fp_cog, fp_so_j );
                sv_add( 3, fp_so_j, rao_j, fp_so_j );

                // Calculate movement of the field point due to the motion of the body
                if ( mode_wl )
                {
                    // Create direction vector
                    tv[0]   = ( panel_geom_->x_wl[1] - panel_geom_->x_wl[0] ) / panel_geom_->len_wl;
                    tv[1]   = ( panel_geom_->y_wl[1] - panel_geom_->y_wl[0] ) / panel_geom_->len_wl;

                    copy_vector( 3, vel_j, v0 );
                    conj_vector( 3, v0, v0 );
                    cross( v0, fp_so_i, v1 );

                    
                    qb_rhs[index]   -= 0.25 * (
                                                    v1[0] * tv[0] + v1[1] * tv[1] + v1[2] * 0.0
                                                ) * G[ fpi ];

                    copy_vector( 3, fp_so_i, v0 );
                    conj_vector( 3, v0, v0 );
                    cross( vel_j, v0, v1 );

                    
                    qb_rhs[index]   -= 0.25 * (
                                                    v1[0] * tv[0] + v1[1] * tv[1] + v1[2] * 0.0
                                                ) * G[ fpi ];

                }
                else
                {
                    // Calculate wave exciting contribution
                    wave_dx         = wave_potential_so_space_dx( 
                                                                        field_point[0],
                                                                        field_point[1],
                                                                        field_point[2],
                                                                        wd,
                                                                        qtf_type
                                                                    );
    
                    wave_dy         = wave_potential_so_space_dy( 
                                                                        field_point[0],
                                                                        field_point[1],
                                                                        field_point[2],
                                                                        wd,
                                                                        qtf_type
                                                                    );
    
                    wave_dz         = wave_potential_so_space_dz( 
                                                                        field_point[0],
                                                                        field_point[1],
                                                                        field_point[2],
                                                                        wd,
                                                                        qtf_type
                                                                    );
    
                    qb_rhs[index]   += -(
                                                wave_dx * nvec[0]
                                                +
                                                wave_dy * nvec[1]
                                                +
                                                wave_dz * nvec[2]
                                            ) * G[ fpi ];
    
                    // Calculate second order rotation term
                    _calculate_so_rotation_mat( 
                                                    &rao_i[3],
                                                    &rao_j[3],
                                                    qtf_type,
                                                    so_rotmat
                                                );
    
                    qb_rhs[index]   +=  (
                                            nvec[0] * ( so_rotmat[0] * fp_cog[0] )
                                            +
                                            nvec[1] * ( so_rotmat[3] * fp_cog[0] + so_rotmat[4] * fp_cog[1] )
                                            +
                                            nvec[2] * ( so_rotmat[6] * fp_cog[0] + so_rotmat[7] * fp_cog[1] + so_rotmat[8] * fp_cog[2] )
                                        ) * cuscomplex( 0.0, -wi_wj ) * G[ fpi ];
    
                    // Calculate relative velocity term
                    cross( rao_i, nvec, v0 );
    
                    v1[0] = std::conj( cuscomplex( 0.0, -wj ) * fp_so_j[0] - vel_j[0] );
                    v1[1] = std::conj( cuscomplex( 0.0, -wj ) * fp_so_j[1] - vel_j[1] );
                    v1[2] = std::conj( cuscomplex( 0.0, -wj ) * fp_so_j[2] - vel_j[2] );
    
                    qb_rhs[index]   += 0.25 * dot_product( 3, v0, v1 ) * G[ fpi ];
    
                    cross( rao_j, nvec, v0 );
                    conj_vector( 3, v0, v0 );
    
                    v1[0] = cuscomplex( 0.0, -wi ) * fp_so_i[0] - vel_i[0];
                    v1[1] = cuscomplex( 0.0, -wi ) * fp_so_i[1] - vel_i[1];
                    v1[2] = cuscomplex( 0.0, -wi ) * fp_so_i[2] - vel_i[2];
                    
                    qb_rhs[index]   += 0.25 * dot_product( 3, v0, v1 ) * G[ fpi ];
    
                    // Calculate velocity gradient term 0
                    copy_vector( 3, vel_j, v1 );
                    conj_vector( 3, v1, v1 );
    
                    qb_rhs[index]   -= 0.25 * (
                                                    ( nvec[0] * fp_so_i[0] + nvec[1] * fp_so_i[1] + nvec[2] * fp_so_i[2] )
                                                    *
                                                    ( v1[0] * dG_dx[fpi] + v1[1] * dG_dy[fpi] + v1[2] * dG_dz[fpi] )
                                                );
    
                    copy_vector( 3, fp_so_j, v0 );
                    conj_vector( 3, v0, v0 );
    
                    qb_rhs[index]   -= 0.25 * (
                                                    ( nvec[0] * v0[0] + nvec[1] * v0[1] + nvec[2] * v0[2] )
                                                    *
                                                    ( vel_i[0] * dG_dx[fpi] + vel_i[1] * dG_dy[fpi] + vel_i[2] * dG_dz[fpi] )
                                                );
    
                    // Calculate velocity gradient term 1
                    copy_vector( 3, vel_j, v0 );
                    conj_vector( 3, v0, v0 );
                    copy_vector( 3, &(rao_i[3]), v1 );
    
                    qb_rhs[index]   -= 0.25 * ( 
                                                    nvec[0] * ( v0[1] * v1[2] - v0[2] * v1[1] )
                                                    +
                                                    nvec[1] * ( v0[2] * v1[0] - v0[0] * v1[2] )
                                                    +
                                                    nvec[2] * ( v0[0] * v1[1] - v0[1] * v1[0] )
                                                ) * G[ fpi ];
    
                    copy_vector( 3, vel_i, v0 );
                    copy_vector( 3, &(rao_j[3]), v1 );
                    conj_vector( 3, v1, v1 );
    
                    qb_rhs[index]   -= 0.25 * ( 
                                                    nvec[0] * ( v0[1] * v1[2] - v0[2] * v1[1] )
                                                    +
                                                    nvec[1] * ( v0[2] * v1[0] - v0[0] * v1[2] )
                                                    +
                                                    nvec[2] * ( v0[0] * v1[1] - v0[1] * v1[0] )
                                                ) * G[ fpi ];
    
                    // Calculate velocity gradient term 2
                    dvdn            = std::conj( vel_j[0] * nvec[0] + vel_j[1] * nvec[1] + vel_j[2] * nvec[2] );
                    qb_rhs[index]   += 0.25 * dvdn * ( 
                                                            fp_so_i[0] * dG_dx[fpi] 
                                                            + 
                                                            fp_so_i[1] * dG_dy[fpi] 
                                                            + 
                                                            fp_so_i[2] * dG_dz[fpi] 
                                                        );
    
                    dvdn            = vel_i[0] * nvec[0] + vel_i[1] * nvec[1] + vel_i[2] * nvec[2];
                    qb_rhs[index]   += 0.25 * dvdn * ( 
                                                            std::conj( fp_so_j[0] ) * dG_dx[fpi]
                                                            +
                                                            std::conj( fp_so_j[1] ) * dG_dy[fpi]
                                                            +
                                                            std::conj( fp_so_j[2] ) * dG_dz[fpi]
                                                        );
                }
            }
        }
    }
}


template<typename Config>
inline  void    so_potential_qf_rhs(
                                            QTFTypeE                        qtf_type,
                                            Input*                          input,
                                            std::size_t                     freq_i,
                                            std::size_t                     freq_j,
                                            cuscomplex*                     raos_hist,
                                            std::size_t                     body_id,
                                            cusfloat*                       cog,
                                            PanelData<cuscomplex, Config>*  panel_data,
                                            cuscomplex*                     G,
                                            cuscomplex*                     dG_dx,
                                            cuscomplex*                     dG_dy,
                                            cuscomplex*                     dG_dz,
                                            WaveDispersionSO*               wd,
                                            BoundarySubtypeE                bsubtype,
                                            cuscomplex*                     qf_rhs
                                )
{
    // Get local variables
    std::size_t bodies_np       = static_cast<input->bodies_np>;
    std::size_t dofs_np         = static_cast<input->dofs_np>;
    std::size_t heads_np        = static_cast<input->heads_np>;
    cuscomplex  dvdn            = 0.0;
    cusfloat    field_point[3]  = { 0.0, 0.0, 0.0 };
    cuscomplex  fp_cog[3]       = { 0.0, 0.0, 0.0 };
    cuscomplex  fp_so_i[3]      = { 0.0, 0.0, 0.0 };
    cuscomplex  fp_so_j[3]      = { 0.0, 0.0, 0.0 };
    std::size_t index           = 0;
    cusfloat    nlvec[3]        = { 0.0, 0.0, 0.0 };
    cusfloat    nlvec_mod       = 0.0;
    cusfloat*   nvec            = panel_data->panel_geom->normal_vec;
    std::size_t pd_idx          = 0;
    std::size_t pd_jdx          = 0;
    cuscomplex  phi_i           = 0.0;
    cuscomplex  phi_j           = 0.0;
    cuscomplex* rao_i           = &raos_freq_i[ freq_i * heads_np * bodies_np * dofs_np ];
    cuscomplex* rao_j           = &raos_freq_j[ freq_j * heads_np * bodies_np * dofs_np ];
    std::size_t rao_i_idx       = 0;
    std::size_t rao_j_idx       = 0;
    cuscomplex  so_rotmat[9]    ;
    cusfloat    tv[3]           = { 0.0, 0.0, 0.0 };
    cusfloat    tv_mod          = 0.0;
    cuscomplex  vel_i[3]        = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_j[3]        = { 0.0, 0.0, 0.0 };
    cuscomplex  v0[3]           = { 0.0, 0.0, 0.0 };
    cuscomplex  v1[3]           = { 0.0, 0.0, 0.0 };
    cuscomplex  wave_dx         = 0.0;
    cuscomplex  wave_dy         = 0.0;
    cuscomplex  wave_dz         = 0.0;

    // Define QTF sign based on the type of QTF term
    int qtf_sign                = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? -1 : 1;

    // Get angular frequencies
    cusfloat wi                 = input->ang_freqs[freq_i];
    cusfloat wj                 = input->ang_freqs[freq_j];
    cusfloat wi_wj              = wi + qtf_sign * wj;

    // Clear rhs vector to sum up all the constributions on it
    std::size_t size_rhs        = pow2s( heads_np ) * panel_data->field_points_np

    // Loop over headings to get all the possible combinations
    for ( std::size_t ih=0; ih<heads_np; ih++ )
    {
        for ( std::size_t jh=0; jh<heads_np; jh++ )
        {
            // Get RAO values for the current frequency and heading combination
            rao_i_idx   = freq_i * heads_np * bodies_np * dofs_np + ih * bodies_np * dofs_np + body_id * dofs_np;
            rao_j_idx   = freq_j * heads_np * bodies_np * dofs_np + jh * bodies_np * dofs_np + body_id * dofs_np;
            rao_i       = &raos_freq_i[ rao_i_idx ];
            rao_j       = &raos_freq_j[ rao_j_idx ];

            // Update wave dispersion relation for second order wave propagation
            wd->set_new_data( 
                                wi, 
                                wj, 
                                panel_data->head[ih], 
                                panel_data->head[jh] 
                            );

            for ( std::size_t fpi=0; fpi<panel_data->field_points_np; fpi++ )
            {
                index           = ( ih * heads_np + jh ) * panel_data->field_points_np + fpi;
                // Get current field point coordinates
                field_point[0]  = panel_data->field_points[3*fpi+0];
                field_point[1]  = panel_data->field_points[3*fpi+1];
                field_point[2]  = panel_data->field_points[3*fpi+2];

                // Get current field point potential and velocities for both frequencies combinations
                pd_idx          = freq_i * heads_np * panel_data->field_points_np;
                phi_i           = panel_data->pot_total[pd_idx + fpi];
                vel_i[0]        = panel_data->vel_x_total[pd_idx + fpi];
                vel_i[1]        = panel_data->vel_y_total[pd_idx + fpi];
                vel_i[2]        = panel_data->vel_z_total[pd_idx + fpi];
                
                pd_jdx          = freq_j * heads_np * panel_data->field_points_np;
                phi_j           = panel_data->pot_total[pd_jdx + fpi];
                vel_j[0]        = panel_data->vel_x_total[pd_jdx + fpi];
                vel_j[1]        = panel_data->vel_y_total[pd_jdx + fpi];
                vel_j[2]        = panel_data->vel_z_total[pd_jdx + fpi];

                // Define field point coordinates with respect to body cog
                fp_cog[0]       = field_point[0] - cog[0];
                fp_cog[1]       = field_point[1] - cog[1];
                fp_cog[2]       = field_point[2] - cog[2];

                clear_vector( 3, fp_so_i );
                cross( rao_i, fp_cog, fp_so_i );
                sv_add( 3, fp_so_i, rao_i, fp_so_i );

                clear_vector( 3, fp_so_j );
                cross( rao_j, fp_cog, fp_so_j );
                sv_add( 3, fp_so_j, rao_j, fp_so_j );

                // Calculate movement of the field point due to the motion of the body
                if ( 
                        bsubtype == BoundarySubtypeE::WL
                        ||
                        bsubtype == BoundarySubtypeE::PC
                    )
                {
                    // Create direction vector
                    nlvec[0]    = ( panel_geom_->x_wl[1] + panel_geom_->x_wl[0] ) / 2.0;
                    nlvec[1]    = ( panel_geom_->y_wl[1] + panel_geom_->y_wl[0] ) / 2.0;
                    nlvec_mod   = std::sqrt( pow2s( nlvec[0] ) + pow2s( nlvec[1] ) );
                    nlvec[0]    /= nlvec_mod;
                    nlvec[1]    /= nlvec_mod;

                    if ( bsubtype == BoundarySubtypeE::WL )
                    {
                        nlvec[0]    = -nlvec[0];
                        nlvec[1]    = -nlvec[1];
                    }

                    // Calculate normal vector to the waterline panel outwards of the FS surface lid
                    tv[0]       = ( panel_geom_->y_wl[1] - panel_geom_->y_wl[0] ) / panel_geom_->len_wl; // Changed to be the normal with anti-clockwise rotation
                    tv[1]       = ( panel_geom_->x_wl[1] - panel_geom_->x_wl[0] ) / panel_geom_->len_wl; // Changed to be the normal with anti-clockwise rotation

                    // Check normal orientation
                    if ( tv[0] * nlvec[0] + tv[1] * nlvec[1] < 0.0 )
                    {
                        tv[0]   = -tv[0];
                        tv[1]   = -tv[1];
                    }
                    
                    // Calculate WL term
                    copy_vector( 3, vel_j, v0 );
                    conj_vector( 3, v0, v0 );
                    
                    qf_rhs[index]   += cuscomplex( 0.0, 0.25 * wi ) * (
                                                                            v0[0] * tv[0] + v0[1] * tv[1]
                                                                        ) * phi_i * G[ fpi ];

                    qf_rhs[index]   -= cuscomplex( 0.0, 0.25 * wj ) * (
                                                                            vel_i[0] * tv[0] + vel_i[1] * tv[1]
                                                                        ) * std::conj( phi_j ) * G[ fpi ];

                }
                else
                {
                    // Calculate term 1
                    qf_rhs[index]   += (
                                            cuscomplex( 0.0, 0.25 * wi * pow2s( wj ) ) 
                                            * 
                                            phi_i
                                            *
                                            std::conj( vel_j[2] )
                                        ) * G[ fpi ];

                    qf_rhs[index]   -= (
                                            cuscomplex( 0.0, 0.25 * wj * pow2s( wi ) ) 
                                            * 
                                            std::conj( phi_j )
                                            *
                                            vel_i[2]
                                        ) * G[ fpi ];

                    // Calculate term 2
                    copy_vector( 3, vel_j, v0 );
                    conj_vector( 3, v0, v0 );
                    
                    qf_rhs[index]   -= 0.25 * wi * dot_product( 3, vel_i, v0 ) * G[ fpi ];

                    copy_vector( 3, vel_i, v0 );
                    conj_vector( 3, v0, v0 );
                    
                    qf_rhs[index]   += 0.25 * wj * dot_product( 3, vel_j, v0 ) * G[ fpi ];

                    // Calculate term 3
                    copy_vector( 3, vel_j, v0 );
                    conj_vector( 3, v0, v0 );

                    qf_rhs[index]   -= 0.25 * wi * (
                                                        dG_dx[fpi] * v0[0] + dG_dy[fpi] * v0[1] + dG_dz[fpi] * v0[2]
                                                    ) * phi_i;


                    qf_rhs[index]   += 0.25 * wi * (
                                                        dG_dx[fpi] * vel_i[0] + dG_dy[fpi] * vel_i[1] + dG_dz[fpi] * vel_i[2]
                                                    ) * std::conj( phi_j );

                    // Calculate term 4
                    qf_rhs[index]   -= cuscomplex( 0.0, 0.5 * wi_wj ) * dot_product( 3, vel_i, vel_j ) * G[ fpi ];

                }
            }
        }
    }
}