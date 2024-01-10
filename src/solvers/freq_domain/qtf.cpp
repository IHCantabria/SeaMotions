
// Include general usage libraries
#include <cassert>

// Include local modules
#include "qtf.hpp"

#include "froude_krylov.hpp"
#include "../../math/gauss.hpp"
#include "../../math/euler_transforms.hpp"
#include "../../math/math_interface.hpp"
#include "../../math/math_tools.hpp"
#include "../../math/topology.hpp"
#include "../../waves.hpp"


void    calculate_pinkster(
                                Input*      input,
                                MpiConfig*  mpi_config,
                                MeshGroup*  mesh_gp,
                                cusfloat    ang_freq_i,
                                cusfloat    ang_freq_j,
                                cuscomplex* qtf_values
                            )
{
    // Define constants to compute Pinkster force
    cusfloat    ai                              = input->wave_amplitude;
    cusfloat    aj                              = input->wave_amplitude;
    cusfloat    aij                             = 0.0;
    cusfloat    bij                             = 0.0;
    cusfloat    cij                             = 0.0;
    cusfloat    dk                              = 0.0;
    cusfloat    dkwij                           = 0.0;
    cusfloat    fij                             = 0.0;
    cuscomplex  froude_krylov[input->dofs_np];  clear_vector( input->dofs_np, froude_krylov );
    cusfloat    g                               = input->grav_acc;
    cusfloat    ki                              = 0.0;
    cusfloat    kj                              = 0.0;
    cusfloat    h                               = input->water_depth;
    int         idx0                            = 0;
    cusfloat    wij                             = ang_freq_i - ang_freq_j;

    // Clear input vector in order avoid spurious data 
    // to be added to the QTF matrix
    clear_vector(
                    input->heads_np * input->bodies_np * input->dofs_np,
                    qtf_values
                );
    
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
                // Define chunk index
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;

            // Calculate ith, jth wave numbers
            ki = w2k( ang_freq_i, input->water_depth, input->grav_acc );
            kj = w2k( ang_freq_j, input->water_depth, input->grav_acc );

            // Calculate wave numbers difference and the equivalent angular frequency
            dk  = ki - kj;
            dkwij = k2w( dk, input->water_depth, input->grav_acc );

            // Calculate equivalent Froude-Krylov force for frequency difference
            calculate_froude_krylov(
                                        input,
                                        mpi_config,
                                        mesh_gp,
                                        dkwij,
                                        froude_krylov
                                    );
            
            // Calculate correction coefficients
            cij = 2 * ki * kj * wij * ( 1 + std::tanh( ki * h ) * std::tanh( kj * h ) ) / ang_freq_i / ang_freq_j;
            bij = (
                        pow2s( ki ) / ang_freq_i / pow2s( std::cosh( ki * h ) )
                        -
                        pow2s( kj ) / ang_freq_j / pow2s( std::cosh( kj * h ) )
                    );
            aij = (
                        ( bij + cij )
                        /
                        ( pow2s( wij ) - dk * g * std::tanh( dk * h ) )
                    ) * 0.5 * g;
            fij = ai * aj * aij * wij / g;

            // Calcuate second order force contribution using Pinkster approximation
            for ( int k=0; k<input->dofs_np; k++ )
            {
                qtf_values[idx0+k] = fij * froude_krylov[k];
            }
        }
    }
}


cuscomplex  calculate_qtf_diff_term(
                                        cuscomplex c0,
                                        cuscomplex c1
                                    )
{
    return 0.5 * ( 
                    c0 * std::conj( c1 )
                    +
                    std::conj( c0 ) * c1
                );
}


void        calculate_qtf_terms_force(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            int             qtf_type,
                                            cuscomplex*     mdrift_rel_we_i,
                                            cuscomplex*     mdrift_rel_we_j,
                                            cuscomplex*     raos_i,
                                            cuscomplex*     raos_j,
                                            cuscomplex*     vel_x_i,
                                            cuscomplex*     vel_y_i,
                                            cuscomplex*     vel_z_i,
                                            cuscomplex*     vel_x_j,
                                            cuscomplex*     vel_y_j,
                                            cuscomplex*     vel_z_j,
                                            cusfloat        ang_freq_j,
                                            cuscomplex*     qtf_values,
                                            cuscomplex*     qtf_wl,
                                            cuscomplex*     qtf_bern,
                                            cuscomplex*     qtf_acc,
                                            cuscomplex*     qtf_mom,
                                            MLGCmpx*        pot_gp,
                                            MLGCmpx*        vel_gp,
                                            bool            is_multi_head
                                    )
{
    // Asssert if qtf_type is on the range
    bool    assert_test = ( qtf_type == 0 ) | ( qtf_type == 1 );
    assert( assert_test && "qtf_type variable with values differnt to 0 and 1" );

    // Define aux variables to be used along the function
    GaussPoints gp( input->gauss_order );
    int         idx0        = 0;
    int         idx1_i      = 0;
    int         idx1_j      = 0;
    int         idx2        = 0;
    int         ih2_end     = 0;
    int         ih2_start   = 0;
    cuscomplex  int_mod     = 0.0;
    PanelGeom*  panel_k     = nullptr;
    int         ngp         = input->gauss_order;
    int         ngpf_1d     = input->gauss_np_factor_1d( );
    int         ngpf        = input->gauss_np_factor_2d( );
    cuscomplex  val_mod     = cuscomplex( 0.0, 0.0 );

    // Set second heading loop bounds according with the multi-heading option
    if ( is_multi_head )
    {
        ih2_start   = 0;
        ih2_end     = input->heads_np;
    }
    else
    {
        ih2_start   = 0;
        ih2_end     = 1;
    }

    // Clear QTF input vector to ensure that previous data will not be storaged
    // erroneously
    int heads_factor_np = input->heads_np;
    if ( is_multi_head )
    {
        heads_factor_np *= input->heads_np;
    }

    clear_vector( 
                    heads_factor_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_values
                );
    clear_vector( 
                    heads_factor_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_wl
                );
    clear_vector( 
                    heads_factor_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_bern
                );
    clear_vector( 
                    heads_factor_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_acc
                );
    clear_vector( 
                    heads_factor_np * mesh_gp->meshes_np * input->dofs_np,
                    qtf_mom
                );

    // Calculate second order force due to relative wave
    // elevation at the WL
    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                if ( is_multi_head )
                {
                    idx0 = (
                                ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
                                +
                                ih2 * ( input->dofs_np * input->bodies_np )
                                +
                                j * input->dofs_np
                            );
                }
                else
                {
                    idx0 = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                }

                for ( int k=pot_gp->field_points_cnp[j]/ngpf_1d; k<pot_gp->field_points_cnp[j+1]/ngpf_1d; k++ )
                {
                    panel_k = mesh_gp->panels_wl[k];
                    int_mod = cuscomplex( 0.0, 0.0 );

                    for ( int gpi=0; gpi<input->gauss_order; gpi++ )
                    {
                        if ( is_multi_head )
                        {
                            idx1_i = ih1 *  pot_gp->field_points_np + k * ngp + gpi;
                            idx1_j = ih2 *  pot_gp->field_points_np + k * ngp + gpi;
                        }
                        else
                        {
                            idx1_i = ih1 *  pot_gp->field_points_np + k * ngp + gpi;
                            idx1_j = ih1 *  pot_gp->field_points_np + k * ngp + gpi;
                        }

                        if ( qtf_type == 0 )
                        {
                            val_mod = calculate_qtf_diff_term( mdrift_rel_we_i[idx1_i], mdrift_rel_we_j[idx1_j] );
                        }
                        else
                        {
                            val_mod  = mdrift_rel_we_i[idx1_i] * mdrift_rel_we_j[idx1_j];
                        }
                        int_mod += val_mod * gp.weights[gpi] * panel_k->len_wl / 2.0;
                    }


                    for ( int r=0; r<input->dofs_np; r++ )
                    {
                        qtf_values[idx0+r] += - 0.25 * input->grav_acc * input->water_density * int_mod * panel_k->normal_vec[r];
                        qtf_wl[idx0+r]     += - 0.25 * input->grav_acc * input->water_density * int_mod * panel_k->normal_vec[r];
                    }

                }
            }
        }
    }

    // Calculate second order force due to the bernouilly contribution
    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                if ( is_multi_head )
                {
                    idx0 = (
                                ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
                                +
                                ih2 * ( input->dofs_np * input->bodies_np )
                                + 
                                j * input->dofs_np
                            );

                }
                else
                {
                    idx0 = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                }

                for ( int k=vel_gp->field_points_cnp[j]/ngpf; k<vel_gp->field_points_cnp[j+1]/ngpf; k++ )
                {
                    int_mod             = cuscomplex( 0.0, 0.0 );
                    panel_k             = mesh_gp->panels[k];

                    for ( int gpi=0; gpi<ngp; gpi++ )
                    {
                        for ( int gpj=0; gpj<ngp; gpj++ )
                        {
                            if ( is_multi_head )
                            {
                                idx1_i    = ih1 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                                idx1_j    = ih2 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                            }
                            else
                            {
                                idx1_i    = ih1 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                                idx1_j    = ih1 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                            }

                            if ( qtf_type == 0 )
                            {
                                val_mod     =   (
                                                    calculate_qtf_diff_term( vel_x_i[idx1_i], vel_x_j[idx1_j] )
                                                    +
                                                    calculate_qtf_diff_term( vel_y_i[idx1_i], vel_y_j[idx1_j] )
                                                    +
                                                    calculate_qtf_diff_term( vel_z_i[idx1_i], vel_z_j[idx1_j] )
                                                );
                            }
                            else
                            {
                                val_mod     =   (
                                                    vel_x_i[idx1_i] * vel_x_j[idx1_j]
                                                    +
                                                    vel_y_i[idx1_i] * vel_y_j[idx1_j]
                                                    +
                                                    vel_z_i[idx1_i] * vel_z_j[idx1_j]
                                                );
                            }
                            int_mod     +=  val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                            panel_k->num_nodes,
                                                                                                            panel_k->xl,
                                                                                                            panel_k->yl,
                                                                                                            gp.roots[gpi],
                                                                                                            gp.roots[gpj]
                                                                                                        );
                        }

                    }
                    
                    for ( int r=0; r<input->dofs_np; r++ )
                    {
                        qtf_values[idx0+r] += 0.25 * input->water_density * int_mod * panel_k->normal_vec[r];
                        qtf_bern[idx0+r]   += 0.25 * input->water_density * int_mod * panel_k->normal_vec[r];
                    }
                }
            }
        }
    }

    // Calculate second order force due to acceleration term
    cusfloat    cog_to_fp[3];   clear_vector( 3, cog_to_fp );
    cuscomplex  cog_to_fp_c[3]; clear_vector( 3, cog_to_fp );
    cuscomplex  point_disp[3];  clear_vector( 3, point_disp );
    cuscomplex  rao_trans[3];   clear_vector( 3, rao_trans );
    cuscomplex  rao_rot[3];     clear_vector( 3, rao_rot );
    cuscomplex  vel_x_acc;
    cuscomplex  vel_y_acc;
    cuscomplex  vel_z_acc;

    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                // Get index to locate RAO data
                if ( is_multi_head )
                {
                    idx0 = ( 
                                ih1 * ( input->dofs_np * input->bodies_np * input->heads_np )
                                +
                                ih2 * ( input->dofs_np * input->bodies_np )
                                +
                                j * input->dofs_np
                            );
                }
                else
                {
                    idx0 = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                }

                // Get RAO values for the current body
                idx1_i = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                for ( int r=0; r<3; r++ )
                {
                    rao_trans[r]    = raos_i[idx1_i+r] * input->wave_amplitude;
                    rao_rot[r]      = raos_i[idx1_i+3+r] * input->wave_amplitude;
                }

                for ( int k= vel_gp->field_points_cnp[j]/ngpf; k<vel_gp->field_points_cnp[j+1]/ngpf; k++ )
                {
                    // Get handle to the current panel for calculations
                    panel_k = mesh_gp->panels[k];
                    int_mod = cuscomplex( 0.0, 0.0 );

                    for ( int gpi=0; gpi<input->gauss_order; gpi++ )
                    {
                        for ( int gpj=0; gpj<input->gauss_order; gpj++ )
                        {
                            // Define field points index
                            if ( is_multi_head )
                            {
                                idx1_j = ih2 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                            }
                            else
                            {
                                idx1_j = ih1 * vel_gp->field_points_np + k * pow2s( ngp ) + gpi * ngp + gpj;
                            }
                            idx2 = k * pow2s( ngp ) + gpi * ngp + gpj;

                            // Define vector from cog to field point
                            sv_sub( 3, &(vel_gp->field_points[3*idx2]), panel_k->body_cog, cog_to_fp );
                            for ( int r=0; r<3; r++ )
                            {
                                cog_to_fp_c[r]  = cuscomplex( cog_to_fp[r], 0.0 );
                            }

                            // Calculate first order displacement of the panel centre
                            clear_vector( 3, point_disp );

                            cross(
                                        rao_rot,
                                        cog_to_fp_c,
                                        point_disp
                                );
                            sv_add(
                                        3,
                                        point_disp,
                                        rao_trans,
                                        point_disp
                                    );

                            // Get velocity pressure term
                            vel_x_acc   = input->water_density * cuscomplex( 0.0, -ang_freq_j ) * vel_x_j[idx1_j];
                            vel_y_acc   = input->water_density * cuscomplex( 0.0, -ang_freq_j ) * vel_y_j[idx1_j];
                            vel_z_acc   = input->water_density * cuscomplex( 0.0, -ang_freq_j ) * vel_z_j[idx1_j];

                            // Calculate point displacement
                            if ( qtf_type == 0 )
                            {
                                val_mod = 0.5 * (
                                                    calculate_qtf_diff_term( point_disp[0], vel_x_acc )
                                                    +
                                                    calculate_qtf_diff_term( point_disp[1], vel_y_acc )
                                                    +
                                                    calculate_qtf_diff_term( point_disp[2], vel_z_acc )
                                                );
                            }
                            else
                            {
                                val_mod = 0.5 * (
                                                    point_disp[0] * vel_x_acc
                                                    +
                                                    point_disp[1] * vel_y_acc
                                                    +
                                                    point_disp[2] * vel_z_acc
                                                );
                            }
                            int_mod     +=  val_mod * gp.weights[gpi] * gp.weights[gpj] * jacobi_det_2d( 
                                                                                                            panel_k->num_nodes,
                                                                                                            panel_k->xl,
                                                                                                            panel_k->yl,
                                                                                                            gp.roots[gpi],
                                                                                                            gp.roots[gpj]
                                                                                                        );
                        }
                    }

                    for ( int r=0; r<input->dofs_np; r++ )
                    {
                        qtf_values[idx0+r] += int_mod * panel_k->normal_vec[r];
                        qtf_acc[idx0+r]    += int_mod * panel_k->normal_vec[r];
                    }
                }
            }
        }
    }

    // Calculate second order force due to momentum
    cusfloat    ang_2                           = pow2s( ang_freq_j );
    cuscomplex  conj_vec[3];                    clear_vector( 3, conj_vec );
    cuscomplex  hydro_force[input->dofs_np];    clear_vector( input->dofs_np, hydro_force );
    cuscomplex  mom_i[3];                       clear_vector( 3, mom_i );

    for ( int ih1=0; ih1<input->heads_np; ih1++ )
    {
        for ( int ih2=ih2_start; ih2<ih2_end; ih2++ )
        {
            for ( int j=0; j<mesh_gp->meshes_np; j++ )
            {
                // Define chunk index
                idx0    = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;

                if ( is_multi_head )
                {
                    idx1_i  = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                    idx1_j  = ih2 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                }
                else
                {
                    idx1_i  = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                    idx1_j  = ih1 * ( input->dofs_np * input->bodies_np ) + j * input->dofs_np;
                }

                // Calculate total hydrodynamic force
                hydro_force[0] = - input->bodies[j]->mass * raos_j[idx1_i] * ang_2;
                hydro_force[1] = - input->bodies[j]->mass * raos_j[idx1_i+1] * ang_2;
                hydro_force[2] = - input->bodies[j]->mass * raos_j[idx1_i+2] * ang_2;
                hydro_force[3] = - (
                                        input->bodies[j]->inertia[0] * raos_j[idx1_i+3]
                                        +
                                        input->bodies[j]->inertia[1] * raos_j[idx1_i+4]
                                        +
                                        input->bodies[j]->inertia[2] * raos_j[idx1_i+5]
                                    ) * ang_2;
                hydro_force[4] = - (
                                        input->bodies[j]->inertia[1] * raos_j[idx1_i+3]
                                        +
                                        input->bodies[j]->inertia[3] * raos_j[idx1_i+4]
                                        +
                                        input->bodies[j]->inertia[4] * raos_j[idx1_i+5]
                                    ) * ang_2;
                hydro_force[5] = - (
                                        input->bodies[j]->inertia[2] * raos_j[idx1_i+3]
                                        +
                                        input->bodies[j]->inertia[4] * raos_j[idx1_i+4]
                                        +
                                        input->bodies[j]->inertia[5] * raos_j[idx1_i+5]
                                    ) * ang_2;

                // Add moment due to translational forces
                if ( qtf_type == 0 )
                {
                    cuscomplex      scale_f( 0.25, 0.0 );

                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  &(raos_i[idx1_j+3]),    conj_vec                            );
                    cross(          conj_vec,           hydro_force,            mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  hydro_force,            conj_vec                            );
                    cross(          &(raos_i[idx1_j+3]),conj_vec,               mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

                    // Add moment due to rotational force
                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  &(raos_i[idx1_j+3]),    conj_vec                            );
                    cross(          conj_vec,           &(hydro_force[3]),      mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );

                    clear_vector(   3,                  mom_i                                                       );
                    conj_vector(    3,                  &(hydro_force[3]),      conj_vec                            );
                    cross(          &(raos_i[idx1_j+3]),conj_vec,               mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );
                }
                else
                {
                    cuscomplex      scale_f( 0.5, 0.0 );

                    clear_vector(   3,                  mom_i                                                       );
                    cross(          &(raos_i[idx1_j+3]),hydro_force,            mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0]),    mom_i,      &(qtf_values[idx0])     );
                    sv_add(         3,                  &(qtf_mom[idx0]),       mom_i,      &(qtf_mom[idx0])        );

                    // Add moment due to rotational force
                    clear_vector(   3,                  mom_i                                                       );
                    cross(          &(raos_i[idx1_j+3]),&(hydro_force[3]),      mom_i                               );
                    svs_mult(       3,                  mom_i,                  scale_f,    mom_i                   );
                    sv_add(         3,                  &(qtf_values[idx0+3]),  mom_i,      &(qtf_values[idx0+3])   );
                    sv_add(         3,                  &(qtf_mom[idx0+3]),     mom_i,      &(qtf_mom[idx0+3])      );
                }
                
            }
        }
    }
}


void        qtf_distribute_matrix_data(
                                            Input*      input,
                                            int         freq_idx,
                                            int         freq_jdx,
                                            cuscomplex* local_mat,
                                            cuscomplex* global_mat,
                                            int         mode
                                        )
{
    // Declare local variables
    int idx0 = 0;
    int idx1 = 1;

    // Loop over headings to distribute data
    for ( int ih=0; ih<input->heads_np; ih++ )
    {
        for ( int ib=0; ib<input->bodies_np; ib++ )
        {
            for ( int id=0; id<input->dofs_np; id++ )
            {
                idx0    = (
                                ih * ( input->heads_np * input->bodies_np * pow2s( input->angfreqs_np ) * input->dofs_np )
                                +
                                ih * ( input->bodies_np * pow2s( input->angfreqs_np ) * input->dofs_np )
                                +
                                ib * ( pow2s( input->angfreqs_np ) * input->dofs_np )
                                +
                                freq_idx *  ( input->angfreqs_np * input->dofs_np )
                                +
                                freq_jdx *  input->dofs_np
                                +
                                id
                            );
                idx1    = (
                                ih * ( input->bodies_np * input->dofs_np )
                                +
                                ib * input->dofs_np
                                +
                                id
                            );
                
                if ( mode == 0 )
                {
                    global_mat[idx0] = local_mat[idx1];
                }
                else if ( mode == 1 )
                {
                    global_mat[idx0] += local_mat[idx1];
                }
                else
                {
                    std::cout << "ERROR\n function: qtf_distribute_matrix_data - mode: " << mode << " - not valid" << std::endl;
                    throw std::exception( );
                }
            }
        }
    }
}