
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
#include "qtf_indirect_method.hpp"

#include "froude_krylov.hpp"
#include "../../green/kochin.hpp"
#include "../../math/integration.hpp"
#include "../../math/math_tools.hpp"
#include "../../math/special_math.hpp"
#include "../../waves/waves_common.hpp"
#include "tools.hpp"


namespace
{
inline constexpr bool _is_diff_qtf( QTFTypeE qtf_type )
{
    return qtf_type == QTFTypeE::QTF_DIFF_CODE;
}


inline constexpr std::size_t _panel_head_index(
                                                std::size_t freq_pos,
                                                std::size_t heads_np,
                                                std::size_t fp_np,
                                                std::size_t head_id,
                                                std::size_t point_id
                                            )
{
    return freq_pos * ( heads_np * fp_np ) + head_id * fp_np + point_id;
}


inline constexpr std::size_t _panel_mode_index(
                                                std::size_t freq_pos,
                                                std::size_t heads_np,
                                                std::size_t dofs_np,
                                                std::size_t fp_np,
                                                std::size_t mode_id,
                                                std::size_t point_id
                                            )
{
    return freq_pos * ( ( heads_np + dofs_np ) * fp_np ) + mode_id * fp_np + point_id;
}


inline constexpr std::size_t _qtf_force_offset(
                                                std::size_t bodies_np,
                                                std::size_t heads_np,
                                                std::size_t dofs_np,
                                                std::size_t ih1,
                                                std::size_t ih2,
                                                std::size_t body_id
                                            )
{
    return (
                ih1 * ( dofs_np * bodies_np * heads_np )
                +
                ih2 * ( dofs_np * bodies_np )
                +
                body_id * dofs_np
            );
}


inline constexpr std::size_t _rao_hist_index(
                                                std::size_t freq_pos,
                                                std::size_t wave_exc_np,
                                                std::size_t bodies_np,
                                                std::size_t dofs_np,
                                                std::size_t head_id,
                                                std::size_t body_id,
                                                std::size_t dof_id
                                            )
{
    return (
                freq_pos * wave_exc_np
                +
                head_id * ( bodies_np * dofs_np )
                +
                body_id * dofs_np
                +
                dof_id
            );
}
}


void        calculate_qtf_indirect_body_term(
                                                    Input*          input,
                                                    std::size_t     freq_pos_i,
                                                    std::size_t     freq_pos_j,
                                                    QTFTypeE        qtf_type,
                                                    cuscomplex*     raos_hist,
                                                    RDDQTF*         body_gp,
                                                    RDDQTF*         wl_gp,
                                                    cuscomplex*     qtf_body_force
                                                )
{
    clear_vector(
                    pow2s( input->heads_np ) * input->bodies_np * input->dofs_np,
                    qtf_body_force
                );

    if ( body_gp == nullptr || wl_gp == nullptr )
    {
        return;
    }

    const bool          is_diff     = _is_diff_qtf( qtf_type );
    const cusfloat      ang_freq_i  = input->angfreqs[freq_pos_i];
    const cusfloat      ang_freq_j  = input->angfreqs[freq_pos_j];
    const std::size_t   bodies_np   = static_cast<std::size_t>( input->bodies_np );
    const std::size_t   dofs_np     = static_cast<std::size_t>( input->dofs_np );
    const std::size_t   heads_np    = static_cast<std::size_t>( input->heads_np );
    const std::size_t   wave_exc_np = heads_np * bodies_np * dofs_np;
    const cusfloat      rho_w       = input->water_density;
    const cusfloat      w_ds        = ( is_diff ) ? ang_freq_i - ang_freq_j : ang_freq_i + ang_freq_j;

    WaveDispersionSO    wdso(
                                input->wave_amplitude,
                                input->wave_amplitude,
                                ang_freq_i,
                                ang_freq_j,
                                input->heads[0],
                                input->heads[0],
                                input->water_depth,
                                input->grav_acc
                            );

    GaussPoints gp( input->gauss_order );
    cuscomplex  ca_i[3]                    = { 0.0, 0.0, 0.0 };
    cuscomplex  ca_j[3]                    = { 0.0, 0.0, 0.0 };
    cuscomplex  cb_i[3]                    = { 0.0, 0.0, 0.0 };
    cuscomplex  field_point_pos_i[3]       = { 0.0, 0.0, 0.0 };
    cuscomplex  field_point_pos_j[3]       = { 0.0, 0.0, 0.0 };
    cuscomplex  field_point_vel_i[3]       = { 0.0, 0.0, 0.0 };
    cuscomplex  field_point_vel_j[3]       = { 0.0, 0.0, 0.0 };
    cuscomplex  fluid_body_vel_raddif[3]   = { 0.0, 0.0, 0.0 };
    cuscomplex  fluid_body_vel_total_i[3]  = { 0.0, 0.0, 0.0 };
    cuscomplex  fluid_body_vel_total_j[3]  = { 0.0, 0.0, 0.0 };
    cuscomplex  normal_vec_c[3]            = { 0.0, 0.0, 0.0 };
    cuscomplex  psi_ds                     = 0.0;
    cuscomplex  rao_rot_i[3]               = { 0.0, 0.0, 0.0 };
    cuscomplex  rao_rot_j[3]               = { 0.0, 0.0, 0.0 };
    cuscomplex  rao_trans_i[3]             = { 0.0, 0.0, 0.0 };
    cuscomplex  rao_trans_j[3]             = { 0.0, 0.0, 0.0 };
    cuscomplex  val_mod                    = 0.0;
    cuscomplex  wave_inc[3]                = { 0.0, 0.0, 0.0 };

    for ( std::size_t ih1 = 0; ih1 < heads_np; ++ih1 )
    {
        for ( std::size_t ih2 = 0; ih2 < heads_np; ++ih2 )
        {
            for ( std::size_t panel_id = 0; panel_id < body_gp->get_size_local( ); ++panel_id )
            {
                auto*       paneld      = &( body_gp->panel_data[panel_id] );
                PanelGeom*  panel_k     = paneld->panel_geom;
                std::size_t body_id      = paneld->body_id;
                std::size_t fp_np        = paneld->field_points_np;
                std::size_t force_offset = _qtf_force_offset( bodies_np, heads_np, dofs_np, ih1, ih2, body_id );
                std::size_t rao_i_base   = _rao_hist_index( freq_pos_i, wave_exc_np, bodies_np, dofs_np, ih1, body_id, 0 );
                std::size_t rao_j_base   = _rao_hist_index( freq_pos_j, wave_exc_np, bodies_np, dofs_np, ih2, body_id, 0 );

                for ( int r = 0; r < 3; ++r )
                {
                    normal_vec_c[r] = cuscomplex( panel_k->normal_vec[r], 0.0 );
                }

                for ( int gpi = 0; gpi < input->gauss_order; ++gpi )
                {
                    for ( int gpj = 0; gpj < input->gauss_order; ++gpj )
                    {
                        std::size_t local_idx = static_cast<std::size_t>( gpi * input->gauss_order + gpj );
                        if ( local_idx >= fp_np )
                        {
                            continue;
                        }

                        std::size_t idx1_i = _panel_head_index( freq_pos_i, heads_np, fp_np, ih1, local_idx );
                        std::size_t idx1_j = _panel_head_index( freq_pos_j, heads_np, fp_np, ih2, local_idx );
                        cusfloat*   field_point = &( paneld->field_points[3 * local_idx] );

                        wave_inc[0] = wave_potential_so_space_dx( field_point[0], field_point[1], field_point[2], &wdso, static_cast<QTFTypeE>( qtf_type ) );
                        wave_inc[1] = wave_potential_so_space_dy( field_point[0], field_point[1], field_point[2], &wdso, static_cast<QTFTypeE>( qtf_type ) );
                        wave_inc[2] = wave_potential_so_space_dz( field_point[0], field_point[1], field_point[2], &wdso, static_cast<QTFTypeE>( qtf_type ) );

                        val_mod = (
                                        wave_inc[0] * panel_k->normal_vec[0]
                                        +
                                        wave_inc[1] * panel_k->normal_vec[1]
                                        +
                                        wave_inc[2] * panel_k->normal_vec[2]
                                  );

                        for ( int r = 0; r < 3; ++r )
                        {
                            rao_rot_i[r]    = raos_hist[rao_i_base + 3 + r] * input->wave_amplitude;
                            rao_rot_j[r]    = raos_hist[rao_j_base + 3 + r] * input->wave_amplitude;
                            rao_trans_i[r]  = raos_hist[rao_i_base + r]     * input->wave_amplitude;
                            rao_trans_j[r]  = raos_hist[rao_j_base + r]     * input->wave_amplitude;
                        }

                        calculate_field_point_vel_rot(
                                                        rao_trans_i,
                                                        rao_rot_i,
                                                        field_point,
                                                        panel_k->body_cog,
                                                        ang_freq_i,
                                                        field_point_vel_i
                                                    );

                        calculate_field_point_vel_rot(
                                                        rao_trans_j,
                                                        rao_rot_j,
                                                        field_point,
                                                        panel_k->body_cog,
                                                        ang_freq_j,
                                                        field_point_vel_j
                                                    );

                        fluid_body_vel_total_i[0] = paneld->vel_x_total[idx1_i];
                        fluid_body_vel_total_i[1] = paneld->vel_y_total[idx1_i];
                        fluid_body_vel_total_i[2] = paneld->vel_z_total[idx1_i];
                        fluid_body_vel_total_j[0] = paneld->vel_x_total[idx1_j];
                        fluid_body_vel_total_j[1] = paneld->vel_y_total[idx1_j];
                        fluid_body_vel_total_j[2] = paneld->vel_z_total[idx1_j];

                        if ( is_diff )
                        {
                            conj_vector( 3, field_point_vel_j, field_point_vel_j );
                            conj_vector( 3, fluid_body_vel_total_j, fluid_body_vel_total_j );
                        }

                        sv_sub( 3, field_point_vel_i, fluid_body_vel_total_i, ca_i );
                        sv_sub( 3, field_point_vel_j, fluid_body_vel_total_j, ca_j );

                        val_mod -= 0.5 * (
                                                sv_dot( 3, ca_i, field_point_vel_i )
                                                +
                                                sv_dot( 3, ca_j, field_point_vel_j )
                                            );

                        cuscomplex int_mod = val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d(
                                                                                                                panel_k->num_nodes,
                                                                                                                panel_k->xl,
                                                                                                                panel_k->yl,
                                                                                                                gp.roots[gpi],
                                                                                                                gp.roots[gpj]
                                                                                                            );

                        for ( std::size_t dof_id = 0; dof_id < dofs_np; ++dof_id )
                        {
                            std::size_t mode_i = _panel_mode_index( freq_pos_i, heads_np, dofs_np, fp_np, dof_id, local_idx );
                            std::size_t mode_j = _panel_mode_index( freq_pos_j, heads_np, dofs_np, fp_np, dof_id, local_idx );
                            psi_ds = ( is_diff ) ? ( paneld->pot_modes[mode_i] - paneld->pot_modes[mode_j] )
                                                 : ( paneld->pot_modes[mode_i] + paneld->pot_modes[mode_j] );

                            qtf_body_force[force_offset + dof_id] += cuscomplex( 0.0, w_ds * rho_w ) * psi_ds * int_mod;
                        }
                    }
                }
            }
        }
    }

    for ( std::size_t ih1 = 0; ih1 < heads_np; ++ih1 )
    {
        for ( std::size_t ih2 = 0; ih2 < heads_np; ++ih2 )
        {
            for ( std::size_t panel_id = 0; panel_id < body_gp->get_size_local( ); ++panel_id )
            {
                auto*       paneld      = &( body_gp->panel_data[panel_id] );
                PanelGeom*  panel_k     = paneld->panel_geom;
                std::size_t body_id      = paneld->body_id;
                std::size_t fp_np        = paneld->field_points_np;
                std::size_t force_offset = _qtf_force_offset( bodies_np, heads_np, dofs_np, ih1, ih2, body_id );
                std::size_t rao_i_base   = _rao_hist_index( freq_pos_i, wave_exc_np, bodies_np, dofs_np, ih1, body_id, 0 );
                std::size_t rao_j_base   = _rao_hist_index( freq_pos_j, wave_exc_np, bodies_np, dofs_np, ih2, body_id, 0 );

                for ( int r = 0; r < 3; ++r )
                {
                    normal_vec_c[r] = cuscomplex( panel_k->normal_vec[r], 0.0 );
                }

                for ( int gpi = 0; gpi < input->gauss_order; ++gpi )
                {
                    for ( int gpj = 0; gpj < input->gauss_order; ++gpj )
                    {
                        std::size_t local_idx = static_cast<std::size_t>( gpi * input->gauss_order + gpj );
                        if ( local_idx >= fp_np )
                        {
                            continue;
                        }

                        std::size_t idx1_i = _panel_head_index( freq_pos_i, heads_np, fp_np, ih1, local_idx );
                        std::size_t idx1_j = _panel_head_index( freq_pos_j, heads_np, fp_np, ih2, local_idx );
                        cusfloat*   field_point = &( paneld->field_points[3 * local_idx] );

                        for ( int r = 0; r < 3; ++r )
                        {
                            rao_rot_i[r]    = raos_hist[rao_i_base + 3 + r] * input->wave_amplitude;
                            rao_rot_j[r]    = raos_hist[rao_j_base + 3 + r] * input->wave_amplitude;
                            rao_trans_i[r]  = raos_hist[rao_i_base + r]     * input->wave_amplitude;
                            rao_trans_j[r]  = raos_hist[rao_j_base + r]     * input->wave_amplitude;
                        }

                        fluid_body_vel_total_i[0] = paneld->vel_x_total[idx1_i];
                        fluid_body_vel_total_i[1] = paneld->vel_y_total[idx1_i];
                        fluid_body_vel_total_i[2] = paneld->vel_z_total[idx1_i];
                        fluid_body_vel_total_j[0] = paneld->vel_x_total[idx1_j];
                        fluid_body_vel_total_j[1] = paneld->vel_y_total[idx1_j];
                        fluid_body_vel_total_j[2] = paneld->vel_z_total[idx1_j];

                        calculate_field_point_rot(
                                                    rao_trans_i,
                                                    rao_rot_i,
                                                    field_point,
                                                    panel_k->body_cog,
                                                    field_point_pos_i
                                                );

                        calculate_field_point_rot(
                                                    rao_trans_j,
                                                    rao_rot_j,
                                                    field_point,
                                                    panel_k->body_cog,
                                                    field_point_pos_j
                                                );

                        if ( is_diff )
                        {
                            conj_vector( 3, field_point_pos_j, field_point_pos_j );
                            conj_vector( 3, fluid_body_vel_total_j, fluid_body_vel_total_j );
                        }

                        for ( std::size_t dof_id = 0; dof_id < dofs_np; ++dof_id )
                        {
                            std::size_t mode_i = _panel_mode_index( freq_pos_i, heads_np, dofs_np, fp_np, dof_id, local_idx );
                            std::size_t mode_j = _panel_mode_index( freq_pos_j, heads_np, dofs_np, fp_np, dof_id, local_idx );

                            fluid_body_vel_raddif[0] = ( is_diff ) ? ( paneld->vel_x_modes[mode_i] - paneld->vel_x_modes[mode_j] )
                                                                   : ( paneld->vel_x_modes[mode_i] + paneld->vel_x_modes[mode_j] );
                            fluid_body_vel_raddif[1] = ( is_diff ) ? ( paneld->vel_y_modes[mode_i] - paneld->vel_y_modes[mode_j] )
                                                                   : ( paneld->vel_y_modes[mode_i] + paneld->vel_y_modes[mode_j] );
                            fluid_body_vel_raddif[2] = ( is_diff ) ? ( paneld->vel_z_modes[mode_i] - paneld->vel_z_modes[mode_j] )
                                                                   : ( paneld->vel_z_modes[mode_i] + paneld->vel_z_modes[mode_j] );

                            val_mod = (
                                            sv_dot( 3, fluid_body_vel_total_j, fluid_body_vel_raddif ) * sv_dot( 3, field_point_pos_i, normal_vec_c )
                                            +
                                            sv_dot( 3, fluid_body_vel_total_i, fluid_body_vel_raddif ) * sv_dot( 3, field_point_pos_j, normal_vec_c )
                                            -
                                            sv_dot( 3, fluid_body_vel_raddif, field_point_pos_i ) * sv_dot( 3, fluid_body_vel_total_j, normal_vec_c )
                                            -
                                            sv_dot( 3, fluid_body_vel_raddif, field_point_pos_j ) * sv_dot( 3, fluid_body_vel_total_i, normal_vec_c )
                                      );

                            cuscomplex int_mod = val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d(
                                                                                                                    panel_k->num_nodes,
                                                                                                                    panel_k->xl,
                                                                                                                    panel_k->yl,
                                                                                                                    gp.roots[gpi],
                                                                                                                    gp.roots[gpj]
                                                                                                                );

                            qtf_body_force[force_offset + dof_id] += cuscomplex( 0.0, w_ds * rho_w / 2.0 ) * int_mod;
                        }
                    }
                }
            }
        }
    }

    for ( std::size_t ih1 = 0; ih1 < heads_np; ++ih1 )
    {
        for ( std::size_t ih2 = 0; ih2 < heads_np; ++ih2 )
        {
            for ( std::size_t panel_id = 0; panel_id < wl_gp->get_size_local( ); ++panel_id )
            {
                auto*       paneld      = &( wl_gp->panel_data[panel_id] );
                PanelGeom*  panel_k     = paneld->panel_geom;
                std::size_t body_id      = paneld->body_id;
                std::size_t fp_np        = paneld->field_points_np;
                std::size_t force_offset = _qtf_force_offset( bodies_np, heads_np, dofs_np, ih1, ih2, body_id );
                std::size_t rao_i_base   = _rao_hist_index( freq_pos_i, wave_exc_np, bodies_np, dofs_np, ih1, body_id, 0 );
                std::size_t rao_j_base   = _rao_hist_index( freq_pos_j, wave_exc_np, bodies_np, dofs_np, ih2, body_id, 0 );

                for ( int r = 0; r < 3; ++r )
                {
                    normal_vec_c[r] = cuscomplex( panel_k->normal_vec[r], 0.0 );
                }

                for ( int gpi = 0; gpi < input->gauss_order; ++gpi )
                {
                    std::size_t local_idx = static_cast<std::size_t>( gpi );
                    if ( local_idx >= fp_np )
                    {
                        continue;
                    }

                    std::size_t idx1_i = _panel_head_index( freq_pos_i, heads_np, fp_np, ih1, local_idx );
                    std::size_t idx1_j = _panel_head_index( freq_pos_j, heads_np, fp_np, ih2, local_idx );
                    cusfloat*   field_point = &( paneld->field_points[3 * local_idx] );

                    for ( int r = 0; r < 3; ++r )
                    {
                        rao_rot_i[r]    = raos_hist[rao_i_base + 3 + r] * input->wave_amplitude;
                        rao_rot_j[r]    = raos_hist[rao_j_base + 3 + r] * input->wave_amplitude;
                        rao_trans_i[r]  = raos_hist[rao_i_base + r]     * input->wave_amplitude;
                        rao_trans_j[r]  = raos_hist[rao_j_base + r]     * input->wave_amplitude;
                    }

                    fluid_body_vel_total_i[0] = paneld->vel_x_total[idx1_i];
                    fluid_body_vel_total_i[1] = paneld->vel_y_total[idx1_i];
                    fluid_body_vel_total_i[2] = paneld->vel_z_total[idx1_i];
                    fluid_body_vel_total_j[0] = paneld->vel_x_total[idx1_j];
                    fluid_body_vel_total_j[1] = paneld->vel_y_total[idx1_j];
                    fluid_body_vel_total_j[2] = paneld->vel_z_total[idx1_j];

                    calculate_field_point_rot(
                                                rao_trans_i,
                                                rao_rot_i,
                                                field_point,
                                                panel_k->body_cog,
                                                field_point_pos_i
                                            );

                    calculate_field_point_rot(
                                                rao_trans_j,
                                                rao_rot_j,
                                                field_point,
                                                panel_k->body_cog,
                                                field_point_pos_j
                                            );

                    if ( is_diff )
                    {
                        conj_vector( 3, field_point_pos_j, field_point_pos_j );
                        conj_vector( 3, fluid_body_vel_total_j, fluid_body_vel_total_j );
                    }

                    cross( field_point_pos_i, fluid_body_vel_total_j, ca_i );
                    cross( field_point_pos_j, fluid_body_vel_total_i, ca_j );

                    for ( std::size_t dof_id = 0; dof_id < dofs_np; ++dof_id )
                    {
                        std::size_t mode_i = _panel_mode_index( freq_pos_i, heads_np, dofs_np, fp_np, dof_id, local_idx );
                        std::size_t mode_j = _panel_mode_index( freq_pos_j, heads_np, dofs_np, fp_np, dof_id, local_idx );
                        psi_ds = ( is_diff ) ? ( paneld->pot_modes[mode_i] - paneld->pot_modes[mode_j] )
                                             : ( paneld->pot_modes[mode_i] + paneld->pot_modes[mode_j] );

                        sv_add( 3, ca_i, ca_j, cb_i );
                        svs_mult( 3, cb_i, psi_ds, cb_i );
                        val_mod = sv_dot( 3, cb_i, normal_vec_c );

                        cuscomplex int_mod = val_mod * gp.weights[gpi] * panel_k->len_wl / 2.0;
                        qtf_body_force[force_offset + dof_id] += cuscomplex( 0.0, -w_ds * rho_w / 2.0 ) * int_mod;
                    }
                }
            }
        }
    }
}


void        calculate_qtf_indirect_fs_near_term(
                                                    Input*          input,
                                                    std::size_t     freq_pos_i,
                                                    std::size_t     freq_pos_j,
                                                    QTFTypeE        qtf_type,
                                                    RDDQTF*         fs_gp,
                                                    cuscomplex*     qtf_fs_force
                                                )
{
    clear_vector(
                    pow2s( input->heads_np ) * input->bodies_np * input->dofs_np,
                    qtf_fs_force
                );

    if ( fs_gp == nullptr )
    {
        return;
    }

    const bool          is_diff     = _is_diff_qtf( qtf_type );
    const cusfloat      ang_freq_i  = input->angfreqs[freq_pos_i];
    const cusfloat      ang_freq_j  = input->angfreqs[freq_pos_j];
    const std::size_t   bodies_np   = static_cast<std::size_t>( input->bodies_np );
    const std::size_t   dofs_np     = static_cast<std::size_t>( input->dofs_np );
    const std::size_t   heads_np    = static_cast<std::size_t>( input->heads_np );
    const cusfloat      grav_acc    = input->grav_acc;
    const cusfloat      rho_w       = input->water_density;
    const cusfloat      sf          = ( is_diff ) ? 1.0 : -1.0;
    const cusfloat      w_ds        = ( is_diff ) ? ang_freq_i - ang_freq_j : ang_freq_i + ang_freq_j;

    WaveDispersionSO    wdso(
                                input->wave_amplitude,
                                input->wave_amplitude,
                                ang_freq_i,
                                ang_freq_j,
                                input->heads[0],
                                input->heads[0],
                                input->water_depth,
                                input->grav_acc
                            );

    GaussPoints gp( input->gauss_order );
    cuscomplex  f0          = cuscomplex( 0.0, w_ds );
    cuscomplex  f1          = -cuscomplex( 0.0, ang_freq_i / 2.0 / grav_acc );
    cuscomplex  f2          = sf * cuscomplex( 0.0, ang_freq_j / 2.0 / grav_acc );
    cuscomplex  pot_fk_i    = 0.0;
    cuscomplex  pot_fk_j    = 0.0;
    cuscomplex  pot_pert_i  = 0.0;
    cuscomplex  pot_pert_j  = 0.0;
    cuscomplex  pot_total_i = 0.0;
    cuscomplex  pot_total_j = 0.0;
    cuscomplex  psi_ds      = 0.0;
    cuscomplex  val_mod_0   = 0.0;
    cuscomplex  val_mod_1   = 0.0;
    cuscomplex  vel_fk_i[3]    = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_fk_j[3]    = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_pert_i[3]  = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_pert_j[3]  = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_total_i[3] = { 0.0, 0.0, 0.0 };
    cuscomplex  vel_total_j[3] = { 0.0, 0.0, 0.0 };
    cusfloat    w2_i           = pow2s( ang_freq_i );
    cusfloat    w2_j           = pow2s( ang_freq_j );
    cusfloat    k2_i           = pow2s( wdso.k0 );
    cusfloat    k2_j           = pow2s( wdso.k1 );

    for ( std::size_t ih1 = 0; ih1 < heads_np; ++ih1 )
    {
        for ( std::size_t ih2 = 0; ih2 < heads_np; ++ih2 )
        {
            for ( std::size_t panel_id = 0; panel_id < fs_gp->get_size_local( ); ++panel_id )
            {
                auto*       paneld      = &( fs_gp->panel_data[panel_id] );
                PanelGeom*  panel_k     = paneld->panel_geom;
                std::size_t body_id      = paneld->body_id;
                std::size_t fp_np        = paneld->field_points_np;
                std::size_t force_offset = _qtf_force_offset( bodies_np, heads_np, dofs_np, ih1, ih2, body_id );

                for ( int gpi = 0; gpi < input->gauss_order; ++gpi )
                {
                    for ( int gpj = 0; gpj < input->gauss_order; ++gpj )
                    {
                        std::size_t local_idx = static_cast<std::size_t>( gpi * input->gauss_order + gpj );
                        if ( local_idx >= fp_np )
                        {
                            continue;
                        }

                        std::size_t idx1_i = _panel_head_index( freq_pos_i, heads_np, fp_np, ih1, local_idx );
                        std::size_t idx1_j = _panel_head_index( freq_pos_j, heads_np, fp_np, ih2, local_idx );
                        cusfloat*   field_point = &( paneld->field_points[3 * local_idx] );

                        pot_fk_i        = wave_potential_fo_space( input->wave_amplitude, ang_freq_i, wdso.k0, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih1] );
                        pot_fk_j        = wave_potential_fo_space( input->wave_amplitude, ang_freq_j, wdso.k1, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih2] );
                        pot_total_i     = paneld->pot_total[idx1_i];
                        pot_total_j     = paneld->pot_total[idx1_j];
                        pot_pert_i      = pot_total_i - pot_fk_i;
                        pot_pert_j      = pot_total_j - pot_fk_j;

                        vel_fk_i[0]     = wave_potential_fo_space_dx( input->wave_amplitude, ang_freq_i, wdso.k0, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih1] );
                        vel_fk_i[1]     = wave_potential_fo_space_dy( input->wave_amplitude, ang_freq_i, wdso.k0, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih1] );
                        vel_fk_i[2]     = wave_potential_fo_space_dz( input->wave_amplitude, ang_freq_i, wdso.k0, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih1] );

                        vel_fk_j[0]     = wave_potential_fo_space_dx( input->wave_amplitude, ang_freq_j, wdso.k1, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih2] );
                        vel_fk_j[1]     = wave_potential_fo_space_dy( input->wave_amplitude, ang_freq_j, wdso.k1, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih2] );
                        vel_fk_j[2]     = wave_potential_fo_space_dz( input->wave_amplitude, ang_freq_j, wdso.k1, input->water_depth, grav_acc, field_point[0], field_point[1], field_point[2], input->heads[ih2] );

                        vel_total_i[0]  = paneld->vel_x_total[idx1_i];
                        vel_total_i[1]  = paneld->vel_y_total[idx1_i];
                        vel_total_i[2]  = paneld->vel_z_total[idx1_i];
                        vel_total_j[0]  = paneld->vel_x_total[idx1_j];
                        vel_total_j[1]  = paneld->vel_y_total[idx1_j];
                        vel_total_j[2]  = paneld->vel_z_total[idx1_j];

                        vel_pert_i[0]   = vel_total_i[0] - vel_fk_i[0];
                        vel_pert_i[1]   = vel_total_i[1] - vel_fk_i[1];
                        vel_pert_i[2]   = vel_total_i[2] - vel_fk_i[2];
                        vel_pert_j[0]   = vel_total_j[0] - vel_fk_j[0];
                        vel_pert_j[1]   = vel_total_j[1] - vel_fk_j[1];
                        vel_pert_j[2]   = vel_total_j[2] - vel_fk_j[2];

                        if ( is_diff )
                        {
                            pot_fk_j    = std::conj( pot_fk_j );
                            pot_pert_j  = std::conj( pot_pert_j );
                            pot_total_j = std::conj( pot_total_j );
                            conj_vector( 3, vel_fk_j, vel_fk_j );
                            conj_vector( 3, vel_pert_j, vel_pert_j );
                            conj_vector( 3, vel_total_j, vel_total_j );
                        }

                        val_mod_0 = (
                                        f0 * ( sv_dot( 3, vel_total_i, vel_pert_j ) + sv_dot( 3, vel_pert_i, vel_fk_j ) )
                                        +
                                        f1 * ( pot_total_i * ( -w2_j * vel_pert_j[2] ) + pot_pert_i * ( -w2_j * vel_fk_j[2] ) )
                                        +
                                        f2 * ( pot_total_j * ( -w2_i * vel_pert_i[2] ) + pot_pert_j * ( -w2_i * vel_fk_i[2] ) )
                                    );

                        val_mod_1 = 0.5 * (
                                                -cuscomplex( 0.0, ang_freq_i ) * pot_total_i * k2_j * pot_pert_j
                                                -
                                                cuscomplex( 0.0, ang_freq_i ) * pot_pert_i * k2_j * pot_fk_j
                                                +
                                                sf * cuscomplex( 0.0, ang_freq_j ) * pot_total_j * k2_i * pot_pert_i
                                                +
                                                sf * cuscomplex( 0.0, ang_freq_j ) * pot_pert_j * k2_i * pot_fk_i
                                            );

                        cuscomplex int_mod = ( val_mod_0 + val_mod_1 ) * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d(
                                                                                                                                    panel_k->num_nodes,
                                                                                                                                    panel_k->xl,
                                                                                                                                    panel_k->yl,
                                                                                                                                    gp.roots[gpi],
                                                                                                                                    gp.roots[gpj]
                                                                                                                                );

                        for ( std::size_t dof_id = 0; dof_id < dofs_np; ++dof_id )
                        {
                            std::size_t mode_i = _panel_mode_index( freq_pos_i, heads_np, dofs_np, fp_np, dof_id, local_idx );
                            std::size_t mode_j = _panel_mode_index( freq_pos_j, heads_np, dofs_np, fp_np, dof_id, local_idx );
                            psi_ds = ( is_diff ) ? ( paneld->pot_modes[mode_i] - paneld->pot_modes[mode_j] )
                                                 : ( paneld->pot_modes[mode_i] + paneld->pot_modes[mode_j] );

                            qtf_fs_force[force_offset + dof_id] += cuscomplex( 0.0, w_ds * rho_w / grav_acc ) * psi_ds * int_mod;
                        }
                    }
                }
            }
        }
    }
}


void        calculate_qtf_indirect_fs_far_term(
                                                    Input*          input,
                                                    std::size_t     freq_pos_i,
                                                    std::size_t     freq_pos_j,
                                                    QTFTypeE        qtf_type,
                                                    SimulationData* sim_data,
                                                    cuscomplex*     qtf_fs_force
                                                )
{
    // Define local variables
    int         forces_np           = pow2s( input->heads_np ) * input->bodies_np * input->dofs_np;
    int         forces_head_np      = input->bodies_np * input->dofs_np;
    int         idx0                = 0;
    int         idx1                = 0;
    cuscomplex* idf11               = generate_empty_vector<cuscomplex>( forces_np );
    cuscomplex* idf12               = generate_empty_vector<cuscomplex>( forces_np );
    cuscomplex* idf21               = generate_empty_vector<cuscomplex>( forces_np );
    cuscomplex* idf22               = generate_empty_vector<cuscomplex>( forces_np );
    cuscomplex  r0                  = cuscomplex( 0.0, 0.0 );
    cuscomplex  r1                  = cuscomplex( 0.0, 0.0 );
    cuscomplex* th0_aux_force       = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th0_lm1_force       = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th0_l_force         = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th0_lp1_force       = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th1_aux_force       = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th1_lm1_force       = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th1_l_force         = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* th1_lp1_force       = generate_empty_vector<cuscomplex>( forces_head_np );
    cuscomplex* vec_aux_force       = generate_empty_vector<cuscomplex>( forces_head_np );

    // Get required fields
    cusfloat    ang_freq_i          = input->angfreqs[freq_pos_i];
    cusfloat    ang_freq_j          = input->angfreqs[freq_pos_j];

    cuscomplex* kochin_pert_cos_i   = &( sim_data->qtf_kochin_pert_cos_freqs[freq_pos_i*sim_data->qtf_kochin_heads_np] );
    cuscomplex* kochin_pert_sin_i   = &( sim_data->qtf_kochin_pert_sin_freqs[freq_pos_i*sim_data->qtf_kochin_heads_np] );
    cuscomplex* kochin_pert_cos_j   = &( sim_data->qtf_kochin_pert_cos_freqs[freq_pos_j*sim_data->qtf_kochin_heads_np] );
    cuscomplex* kochin_pert_sin_j   = &( sim_data->qtf_kochin_pert_sin_freqs[freq_pos_j*sim_data->qtf_kochin_heads_np] );

    cuscomplex* kochin_rad_cos_i    = &( sim_data->qtf_kochin_rad_cos_freqs[freq_pos_i*sim_data->qtf_kochin_rad_np] );
    cuscomplex* kochin_rad_sin_i    = &( sim_data->qtf_kochin_rad_sin_freqs[freq_pos_i*sim_data->qtf_kochin_rad_np] );
    cuscomplex* kochin_rad_cos_j    = &( sim_data->qtf_kochin_rad_cos_freqs[freq_pos_j*sim_data->qtf_kochin_rad_np] );
    cuscomplex* kochin_rad_sin_j    = &( sim_data->qtf_kochin_rad_sin_freqs[freq_pos_j*sim_data->qtf_kochin_rad_np] );

    cusfloat    grav_acc            = input->grav_acc;
    cusfloat    radius_fs           = input->bodies[0]->mesh_fs_qtf->get_fs_radius( );
    cusfloat    rho_w               = input->water_density;

    // Define second order wave dispersion object
    WaveDispersionSO*   wdso        = new WaveDispersionSO( 
                                                            input->wave_amplitude,
                                                            input->wave_amplitude,
                                                            ang_freq_i,
                                                            ang_freq_j,
                                                            input->heads[0],
                                                            input->heads[0],
                                                            input->water_depth,
                                                            input->grav_acc
                                                        );

    // Define common scaling factors
    cusfloat    k_ds    = 0.0;
    cusfloat    sf0     = 0.0;
    cusfloat    sf1     = 0.0;
    cusfloat    w_ds    = 0.0;

    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        k_ds    = wdso->k_diff_mod;
        sf0     = 1.0;
        sf1     = - 1.0;
        w_ds    = wdso->w_diff;
    }
    else
    {
        k_ds    = wdso->k_sum_mod;
        sf0     = - 1.0;
        sf1     = 1.0;
        w_ds    = wdso->w_sum;
    }

    cuscomplex  sf2     = cuscomplex( 0.0, 1.0 ) * w_ds * rho_w / grav_acc;
    cuscomplex  sf3     = (
                                cuscomplex( 0.0, 1.0 )
                                *
                                input->wave_amplitude
                                *
                                grav_acc
                                *
                                8.0 
                                * 
                                PI
                                *
                                std::sqrt( k_ds )
                                *
                                wave_vertical_profile_mod_fo( 
                                                                k_ds,
                                                                input->water_depth,
                                                                0.0
                                                            )
                            );

    cuscomplex  kappa_1 = sf0 * cuscomplex( 0.0, -1.0 ) * w_ds * wdso->k0 * wdso->k1;
    cuscomplex  kappa_2 =  (
                                (
                                    cuscomplex( 0.0, 1.0 )
                                    *
                                    w_ds
                                    *
                                    pow2s( ang_freq_i ) / grav_acc
                                    *
                                    pow2s( ang_freq_j ) / grav_acc
                                )
                                +
                                sf0 * (
                                            cuscomplex( 0.0, 1.0 ) * ang_freq_i * ang_freq_j / 2.0
                                            *
                                            (
                                                pow2s( wdso->k0 ) / ( ang_freq_i * pow2s( std::cosh( wdso->k0 * input->water_depth ) ) )
                                                +
                                                sf1 * pow2s( wdso->k1 ) / ( ang_freq_j * pow2s( std::cosh( wdso->k1 * input->water_depth ) ) )
                                            )
                                        )
                            );

    // Calculate second order forces for the different headings combinations
    for ( int ih0=0; ih0<input->heads_np; ih0++ )
    {
        for ( int ih1=0; ih1<input->heads_np; ih1++ )
        {
            // Get index for the storage of the body force
            idx0 = (
                        ih0 * ( input->dofs_np * input->bodies_np * input->heads_np )
                        +
                        ih1 * ( input->dofs_np * input->bodies_np )
                    );

            // Get first theta integrals
            idx1    = ih0 * ( input->bodies_np * input->kochin_np );
            calculate_theta_integral(
                                        input,
                                        input->heads[ih0],
                                        -1,
                                        qtf_type,
                                        0,
                                        &( kochin_pert_cos_j[idx1] ),
                                        &( kochin_pert_sin_j[idx1] ),
                                        kochin_rad_cos_i,
                                        kochin_rad_sin_i,
                                        kochin_rad_cos_j,
                                        kochin_rad_sin_j,
                                        th0_lm1_force
                                    );

            calculate_theta_integral(
                                        input,
                                        input->heads[ih0],
                                        0,
                                        qtf_type,
                                        0,
                                        &( kochin_pert_cos_j[idx1] ),
                                        &( kochin_pert_sin_j[idx1] ),
                                        kochin_rad_cos_i,
                                        kochin_rad_sin_i,
                                        kochin_rad_cos_j,
                                        kochin_rad_sin_j,
                                        th0_l_force
                                    );

            idx1    = ih1 * ( input->bodies_np * input->kochin_np );
            calculate_theta_integral(
                                        input,
                                        input->heads[ih1],
                                        -1,
                                        qtf_type,
                                        1,
                                        &( kochin_pert_cos_i[idx1] ),
                                        &( kochin_pert_sin_i[idx1] ),
                                        kochin_rad_cos_i,
                                        kochin_rad_sin_i,
                                        kochin_rad_cos_j,
                                        kochin_rad_sin_j,
                                        th1_lm1_force
                                    );
            
            calculate_theta_integral(
                                        input,
                                        input->heads[ih1],
                                        0,
                                        qtf_type,
                                        1,
                                        &( kochin_pert_cos_i[idx1] ),
                                        &( kochin_pert_sin_i[idx1] ),
                                        kochin_rad_cos_i,
                                        kochin_rad_sin_i,
                                        kochin_rad_cos_j,
                                        kochin_rad_sin_j,
                                        th1_l_force
                                    );

            // Loop over l order in order to perform the sumation
            for ( int l_order=0; l_order<100; l_order++ )
            {
                // Calculate r1 function
                r0 = calculate_r0_integral( radius_fs, wdso, l_order, qtf_type );

                // Calculate r2 function
                r1 = calculate_r1_integral( radius_fs, wdso, l_order, qtf_type );

                // Calculate theta function
                idx1    = ih0 * ( input->bodies_np * input->kochin_np );
                calculate_theta_integral(
                                            input,
                                            input->heads[ih0],
                                            l_order,
                                            qtf_type,
                                            0,
                                            &( kochin_pert_cos_j[idx1] ),
                                            &( kochin_pert_sin_j[idx1] ),
                                            kochin_rad_cos_i,
                                            kochin_rad_sin_i,
                                            kochin_rad_cos_j,
                                            kochin_rad_sin_j,
                                            th0_lp1_force
                                        );

                idx1    = ih1 * ( input->bodies_np * input->kochin_np );
                calculate_theta_integral(
                                            input,
                                            input->heads[ih1],
                                            l_order,
                                            qtf_type,
                                            1,
                                            &( kochin_pert_cos_i[idx1] ),
                                            &( kochin_pert_sin_i[idx1] ),
                                            kochin_rad_cos_i,
                                            kochin_rad_sin_i,
                                            kochin_rad_cos_j,
                                            kochin_rad_sin_j,
                                            th1_lp1_force
                                        );

                // Calculate theta auxiliar function
                sv_add( forces_head_np, th0_lm1_force, th0_lp1_force, th0_aux_force );
                sv_add( forces_head_np, th1_lm1_force, th1_lp1_force, th1_aux_force );

                // Scale theta integrals and sum to their corresponding integrals
                svs_mult( forces_head_np, th0_aux_force, 0.5 * r0, vec_aux_force );
                sv_add( forces_head_np, vec_aux_force, &( idf11[idx0] ), &( idf11[idx0] ));

                svs_mult( forces_head_np, th1_aux_force, 0.5 * r1, vec_aux_force );
                sv_add( forces_head_np, vec_aux_force, &( idf12[idx0] ), &( idf12[idx0] ));

                svs_mult( forces_head_np, th0_l_force, r0, vec_aux_force );
                sv_add( forces_head_np, vec_aux_force, &( idf21[idx0] ), &( idf21[idx0] ));

                svs_mult( forces_head_np, th1_l_force, r1, vec_aux_force );
                sv_add( forces_head_np, vec_aux_force, &( idf22[idx0] ), &( idf22[idx0] ));

                // Roll theta values
                copy_vector( forces_head_np, th0_l_force, th0_lm1_force );
                copy_vector( forces_head_np, th0_lp1_force, th0_l_force );
                copy_vector( forces_head_np, th1_l_force, th1_lm1_force );
                copy_vector( forces_head_np, th1_lp1_force, th1_l_force );

            }

        }
    }

    // Scale radial and angular integrals
    cuscomplex  sf_11 = (
                            -1.0
                            *
                            sf2
                            *
                            kappa_1
                            * 
                            sf3 
                            *
                            std::sqrt( wdso->k1 ) 
                            *
                            wave_vertical_profile_mod_fo( wdso->k1, input->water_depth, 0.0 )
                            /
                            wdso->w0
                        );
    
    cuscomplex  sf_12 = (
                            sf0
                            *
                            sf2
                            *
                            kappa_1
                            *
                            sf3 
                            *
                            std::sqrt( wdso->k0 ) 
                            *
                            wave_vertical_profile_mod_fo( wdso->k0, input->water_depth, 0.0 )
                            /
                            wdso->w1
                        );

    cuscomplex  sf_21 = (
                            -1.0
                            *
                            sf2
                            *
                            kappa_2
                            *
                            sf3
                            *
                            std::sqrt( wdso->k1 ) 
                            *
                            wave_vertical_profile_mod_fo( wdso->k1, input->water_depth, 0.0 )
                            /
                            wdso->w0
                        );
    
    cuscomplex  sf_22 = (
                            sf0
                            *
                            sf2
                            *
                            kappa_2
                            *
                            sf3
                            *
                            std::sqrt( wdso->k0) 
                            *
                            wave_vertical_profile_mod_fo( wdso->k0, input->water_depth, 0.0 )
                            /
                            wdso->w1
                        );

    svs_mult( forces_np, idf11, sf_11, idf11 );
    svs_mult( forces_np, idf12, sf_12, idf12 );
    svs_mult( forces_np, idf21, sf_21, idf21 );
    svs_mult( forces_np, idf22, sf_22, idf22 );

    // Add contributions from the different integrals to have the 
    // comple far field FS contribution
    clear_vector( forces_np, qtf_fs_force );
    sv_add( forces_np, qtf_fs_force, idf11, qtf_fs_force );
    sv_add( forces_np, qtf_fs_force, idf12, qtf_fs_force );
    sv_add( forces_np, qtf_fs_force, idf21, qtf_fs_force );
    sv_add( forces_np, qtf_fs_force, idf22, qtf_fs_force );

    // Delete local heap memory
    mkl_free( idf11 );
    mkl_free( idf12 );
    mkl_free( idf21 );
    mkl_free( idf22 );
    mkl_free( th0_aux_force );
    mkl_free( th0_lm1_force );
    mkl_free( th0_l_force   );
    mkl_free( th0_lp1_force );
    mkl_free( th1_aux_force );
    mkl_free( th1_lm1_force );
    mkl_free( th1_l_force   );
    mkl_free( th1_lp1_force );
    mkl_free( vec_aux_force );
    delete wdso; 
}


cuscomplex  calculate_r0_integral(
                                                    cusfloat            R,
                                                    WaveDispersionSO*   wdso,
                                                    int                 l_order,
                                                    QTFTypeE            qtf_type
                                )
{
    // Calculate scaling factors
    cuscomplex sf0( 1.0, 0.0 );
    if ( qtf_type == QTFTypeE::QTF_SUM_CODE )
    {
        sf0 = std::exp( cuscomplex( 0.0, PI / 2.0 ) );
    }

    cusfloat    ep_l    = ep_n( l_order );
    cuscomplex  sfi     = std::pow( cuscomplex( 0.0, 1.0 ), l_order );

    // Calculate alpha and beta parameters
    cusfloat    alpha   = 0.0;
    cusfloat    beta    = wdso->k0;

    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        alpha   = wdso->k_diff_mod - wdso->k1;
    }
    else
    {
        alpha   = wdso->k_sum_mod + wdso->k1;
    }

    // Calculate 0 to inf analytical integral
    cuscomplex  int_value_ana   = besseljn_expi_int( 
                                                        alpha,
                                                        beta,
                                                        static_cast<cusfloat>( l_order )
                                                    );

    // Calculate 0 to R numerical integral
    auto        cos_kenel_fcn   =   [ alpha, beta, l_order ]
                                    ( cusfloat x )
                                    {
                                        return besseljn_cos_kernel( alpha, beta, l_order, x );
                                    };

    auto        sin_kenel_fcn   =   [ alpha, beta, l_order ]
                                    ( cusfloat x )
                                    {
                                        return besseljn_sin_kernel( alpha, beta, l_order, x );
                                    };

    cusfloat    cos_value_num   = romberg_quadrature(
                                                        cos_kenel_fcn,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );
    
    cusfloat    sin_value_num   = romberg_quadrature(
                                                        sin_kenel_fcn,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );

    cuscomplex  int_value_num   = cuscomplex( cos_value_num, sin_value_num );

    // Calculate R to infinite integral
    cuscomplex  int_value       = sf0 * ep_l * sfi * ( int_value_ana - int_value_num );

    return int_value;
}


cuscomplex  calculate_r1_integral(
                                                    cusfloat            R,
                                                    WaveDispersionSO*   wdso,
                                                    int                 l_order,
                                                    QTFTypeE            qtf_type
                                )
{
    // Define scaling factors
    cuscomplex  sf0     = std::exp( cuscomplex( 0.0, PI / 2.0 ) );
    cusfloat    ep_l    = ep_n( l_order );
    cuscomplex  sfi     = std::pow( cuscomplex( 0.0, 1.0 ), l_order );

    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        sfi = std::conj( sfi );
    }

    // Calculate alpha and beta
    cusfloat    alpha   = 0.0;
    cusfloat    beta    = wdso->k1;
    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
    {
        alpha = wdso->k_diff_mod + wdso->k0;
    }
    else
    {
        alpha = wdso->k_sum_mod + wdso->k0;
    }

    // Calculate 0 to inf analytical integral
    cuscomplex  int_value_ana   = besseljn_expi_int( 
                                                        alpha,
                                                        beta,
                                                        static_cast<cusfloat>( l_order )
                                                    );

    // Calculate 0 to R numerical integral
    auto        cos_kenel_fcn   =   [ alpha, beta, l_order ]
                                    ( cusfloat x )
                                    {
                                        return besseljn_cos_kernel( alpha, beta, l_order, x );
                                    };

    auto        sin_kenel_fcn   =   [ alpha, beta, l_order ]
                                    ( cusfloat x )
                                    {
                                        return besseljn_sin_kernel( alpha, beta, l_order, x );
                                    };
    
    cusfloat    cos_value_num   = romberg_quadrature(
                                                        cos_kenel_fcn,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );
    
    cusfloat    sin_value_num   = romberg_quadrature(
                                                        sin_kenel_fcn,
                                                        0.0,
                                                        R,
                                                        1e-6
                                                    );

    cuscomplex  int_value_num   = cuscomplex( cos_value_num, sin_value_num );

    // Calculate R to infinite integral
    cuscomplex  int_value       = sf0 * ep_l * sfi * ( int_value_ana - int_value_num );

    return int_value;
}


void  		calculate_theta_integral(
                                                    Input*      input,
                                                    cusfloat    beta,
                                                    int         l_order,
                                                    QTFTypeE    qtf_type,
                                                    int         theta_type,
                                                    cuscomplex* kochin_cos_pert_j,
                                                    cuscomplex* kochin_sin_pert_j,
                                                    cuscomplex* kochin_cos_rad_i,
                                                    cuscomplex* kochin_sin_rad_i,
                                                    cuscomplex* kochin_cos_rad_j,
                                                    cuscomplex* kochin_sin_rad_j,
                                                    cuscomplex* body_force
                                    )
{
    // Clear incoming vector in order to avoid taking into account spurious data
    clear_vector( 
                    input->bodies_np * input->dofs_np,
                    body_force
                );

    // Define local variables
    cuscomplex  cpm     = cuscomplex( 0.0, 0.0 );
    cuscomplex  spm     = cuscomplex( 0.0, 0.0 );
    cuscomplex  crn     = cuscomplex( 0.0, 0.0 );
    cuscomplex  srn     = cuscomplex( 0.0, 0.0 );
    int         idx0    = 0;
    int         idx1    = 0;
    cusfloat    t0_int  = 0.0;
    cusfloat    t1_int  = 0.0;
    cusfloat    t2_int  = 0.0;
    cusfloat    t3_int  = 0.0;

    // Loop over values to calculate the resulting forces
    for ( int ib=0; ib<input->bodies_np; ib++ )
    {
        for ( int id=0; id<input->dofs_np; id++ )
        {
            // Define index to stoage the data
            idx0    = ( 
                            ib * input->dofs_np
                            +
                            id
                        );

            // Loop over kochin coefficients to calculate the force
            for ( int m=0; m<input->kochin_np; m++ )
            {
                // Get current index for the mth perturbation coefficient
                idx1 = (
                            ib * input->kochin_np
                            +
                            m
                        );
                
                // Get perturbation series coefficients
                if ( 
                        ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) 
                        && 
                        ( theta_type == 0 )
                    )
                {
                    cpm = std::conj( kochin_cos_pert_j[idx1] );
                    spm = std::conj( kochin_sin_pert_j[idx1] );
                }
                else
                {
                    cpm = kochin_cos_pert_j[idx1];
                    spm = kochin_sin_pert_j[idx1];
                }

                
                for ( int n=0; n<input->kochin_np; n++ )
                {
                    // Get current index for the nth perturbation coefficient
                    idx1 = (
                                id *  ( input->bodies_np * input->kochin_np )
                                +
                                ib * input->kochin_np
                                +
                                m
                            );

                    // Get radiation coefficient
                    if ( qtf_type == QTFTypeE::QTF_DIFF_CODE )
                    {
                        crn = ( kochin_cos_rad_i[idx1] - kochin_cos_rad_j[idx1] );
                        srn = ( kochin_sin_rad_i[idx1] - kochin_sin_rad_j[idx1] );
                    }
                    else
                    {
                        crn = ( kochin_cos_rad_i[idx1] + kochin_cos_rad_j[idx1] );
                        srn = ( kochin_sin_rad_i[idx1] + kochin_sin_rad_j[idx1] );
                    }

                    // Calculate angular integrals
                    t0_int  = calculate_kochin_cosexp_t0( beta, l_order, m, n );
                    t1_int  = calculate_kochin_cosexp_t1( beta, l_order, m, n );
                    t2_int  = calculate_kochin_cosexp_t2( beta, l_order, m, n );
                    t3_int  = calculate_kochin_cosexp_t3( beta, l_order, m, n );

                    // Calculate body force
                    body_force[idx0] += (
                                            cpm * crn * t0_int
                                            +
                                            cpm * srn * t1_int
                                            +
                                            spm * crn * t2_int
                                            +
                                            spm * srn * t3_int
                                        );

                }
            }
        }
    }
}


void        calculate_secord_force_indirect(
                                                    Input*          input,
                                                    MpiConfig*      mpi_config,
                                                    MeshGroup*      mesh_gp,
                                                    std::size_t     freq_pos_i,
                                                    std::size_t     freq_pos_j,
                                                    QTFTypeE        qtf_type,
                                                    RDDQTF*         body_gp,
                                                    RDDQTF*         fs_gp,
                                                    RDDQTF*         wl_gp,
                                                    SimulationData* sim_data
                                            )
{
    const cusfloat ang_freq_i = input->angfreqs[freq_pos_i];
    const cusfloat ang_freq_j = input->angfreqs[freq_pos_j];
    const bool     should_compute =
                        ( ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) && ( freq_pos_i != freq_pos_j ) )
                        ||
                        ( qtf_type == QTFTypeE::QTF_SUM_CODE );

    cuscomplex* qtf_fk_force = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? sim_data->qtf_diff_froude_krylov_fo_p0 : sim_data->qtf_sum_froude_krylov_fo_p0;
    cuscomplex* qtf_body_force = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? sim_data->qtf_diff_body_force_p0 : sim_data->qtf_sum_body_force_p0;
    cuscomplex* qtf_fs_near_force = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? sim_data->qtf_diff_fs_near_field_p0 : sim_data->qtf_sum_fs_near_field_p0;
    cuscomplex* qtf_fs_far_force = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? sim_data->qtf_diff_fs_far_field_p0 : sim_data->qtf_sum_fs_far_field_p0;
    cuscomplex* qtf_second_force = ( qtf_type == QTFTypeE::QTF_DIFF_CODE ) ? sim_data->qtf_diff_secord_force : sim_data->qtf_sum_secord_force;

    if ( mpi_config->is_root( ) )
    {
        clear_vector( sim_data->qtf_np, qtf_fk_force );
        clear_vector( sim_data->qtf_np, qtf_body_force );
        clear_vector( sim_data->qtf_np, qtf_fs_near_force );
        clear_vector( sim_data->qtf_np, qtf_fs_far_force );
        clear_vector( sim_data->qtf_np, qtf_second_force );
    }

    std::vector<cuscomplex> body_force_local( sim_data->qtf_np, cuscomplex( 0.0, 0.0 ) );
    std::vector<cuscomplex> fs_near_force_local( sim_data->qtf_np, cuscomplex( 0.0, 0.0 ) );

    if ( should_compute )
    {
        calculate_qtf_indirect_body_term(
                                            input,
                                            freq_pos_i,
                                            freq_pos_j,
                                            qtf_type,
                                            sim_data->raos_hist,
                                            body_gp,
                                            wl_gp,
                                            body_force_local.data( )
                                        );

        calculate_qtf_indirect_fs_near_term(
                                                input,
                                                freq_pos_i,
                                                freq_pos_j,
                                                qtf_type,
                                                fs_gp,
                                                fs_near_force_local.data( )
                                            );
    }

    cuscomplex* body_reduce_dst = ( mpi_config->is_root( ) ) ? qtf_body_force : body_force_local.data( );
    cuscomplex* fs_reduce_dst   = ( mpi_config->is_root( ) ) ? qtf_fs_near_force : fs_near_force_local.data( );

    MPI_Reduce(
                    body_force_local.data( ),
                    body_reduce_dst,
                    sim_data->qtf_np,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    MPI_Reduce(
                    fs_near_force_local.data( ),
                    fs_reduce_dst,
                    sim_data->qtf_np,
                    mpi_cuscomplex,
                    MPI_SUM,
                    mpi_config->proc_root,
                    MPI_COMM_WORLD
                );

    if ( mpi_config->is_root( ) )
    {
        if ( should_compute )
        {
            calculate_froude_krylov_so(
                                            input,
                                            mesh_gp,
                                            ang_freq_i,
                                            ang_freq_j,
                                            qtf_type,
                                            qtf_fk_force
                                        );

            calculate_qtf_indirect_fs_far_term(
                                                    input,
                                                    freq_pos_i,
                                                    freq_pos_j,
                                                    qtf_type,
                                                    sim_data,
                                                    qtf_fs_far_force
                                                );
        }

        sv_add( sim_data->qtf_np, qtf_second_force, qtf_fk_force, qtf_second_force );
        sv_add( sim_data->qtf_np, qtf_second_force, qtf_body_force, qtf_second_force );
        sv_add( sim_data->qtf_np, qtf_second_force, qtf_fs_near_force, qtf_second_force );
        sv_add( sim_data->qtf_np, qtf_second_force, qtf_fs_far_force, qtf_second_force );
    }
}