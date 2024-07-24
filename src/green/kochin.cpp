
// Include local modules
#include "kochin.hpp"

#include "../math/integration.hpp"
#include "../waves/wave_dispersion_base_fo.hpp"


void        calculate_kochin_coefficients(
                                            Input*              input,
                                            MeshGroup*          mesh_gp,
                                            KochinInterface*    kochin,
                                            cuscomplex*         sources,
                                            cuscomplex*         cos_coeff,
                                            cuscomplex*         sin_coeff
                                        )
{
    // Define local variables
    int         idx0        = 0;
    SourceNode* source_k    = nullptr;

    // Create lambda interface function to perform the adaptive quadrature
    auto interf_fcn =   [ kochin ]
                        ( cusfloat xi, cusfloat eta, cusfloat x, cusfloat y, cusfloat z ) -> cuscomplex
                        {
                            return (*kochin)( xi, eta, x, y, z );
                        };

    // Loop over Kochin expansion series coefficients
    for ( int i=0; i<input->kochin_np; i++ )
    {
        // Set L order to calculate the appropiate coefficient expansion
        kochin->set_l_order( i );
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            // Get index to storage the integral
            idx0    = (
                            j * input->kochin_np
                            +
                            i
                        );

            // Clear coefficients vector memory to not have 
            // spurious data to take into account
            cos_coeff[idx0]  = 0.0;

            for ( int k=mesh_gp->source_nodes_cnp[j]; k<mesh_gp->source_nodes_cnp[j+1]; k++ )
            {
                // Get current panel
                source_k  = mesh_gp->source_nodes[k];

                // Set appropiate source for the integration
                kochin->set_source( source_k );

                // Integrate for cosine coefficient
                kochin->set_is_cos( true );

                cos_coeff[idx0] += sources[k] * adaptive_quadrature_panel(
                                                                            source_k->panel,
                                                                            interf_fcn,
                                                                            input->pot_abs_err,
                                                                            input->pot_rel_err,
                                                                            input->is_block_adaption,
                                                                            false,
                                                                            input->gauss_order
                                                                        );

                // Integrate for sine coefficient
                kochin->set_is_cos( false );

                sin_coeff[idx0] += sources[k] * adaptive_quadrature_panel(
                                                                            source_k->panel,
                                                                            interf_fcn,
                                                                            input->pot_abs_err,
                                                                            input->pot_rel_err,
                                                                            input->is_block_adaption,
                                                                            false,
                                                                            input->gauss_order
                                                                        );
            }
        }
    }
}


cusfloat    calculate_kochin_cosexp_t0(
                                            cusfloat    beta,
                                            int         ln,
                                            int         m,
                                            int         n
                                        )
{
    // Calculate cosine term
    cusfloat    t_cos   = std::cos( ln * beta ) * cos3_int_0_2PI( ln, m, n );

    // Calculate sine term
    cusfloat    t_sin   = std::sin( ln * beta ) * cos2sin_int_0_2PI( m, n, ln );

    return ( t_cos + t_sin );
}


cusfloat    calculate_kochin_cosexp_t1(
                                            cusfloat    beta,
                                            cusfloat    ln,
                                            cusfloat    m,
                                            cusfloat    n
                                        )
{
    // Calculate cosine term
    cusfloat    t_cos   = std::cos( ln * beta ) * cos2sin_int_0_2PI( ln, m, n );

    // Calculate sine term
    cusfloat    t_sin   = std::sin( ln * beta ) * cossin2_int_0_2PI( m, ln, n );

    return ( t_cos + t_sin );
}


cusfloat    calculate_kochin_cosexp_t2(
                                            cusfloat    beta,
                                            cusfloat    ln,
                                            cusfloat    m,
                                            cusfloat    n
                                        )
{
    // Calculate cosine term
    cusfloat    t_cos   = std::cos( ln * beta ) * cos2sin_int_0_2PI( ln, n, m );

    // Calculate sine term
    cusfloat    t_sin   = std::sin( ln * beta ) * cossin2_int_0_2PI( n, ln, m );

    return ( t_cos + t_sin );
}


cusfloat    calculate_kochin_cosexp_t3(
                                            cusfloat    beta,
                                            cusfloat    ln,
                                            cusfloat    m,
                                            cusfloat    n
                                        )
{
    // Calculate cosine term
    cusfloat    t_cos   = std::cos( ln * beta ) * cossin2_int_0_2PI( ln, m, n );

    // Calculate sine term
    cusfloat    t_sin   = std::sin( ln * beta ) * sin3_int_0_2PI( ln, m, n);

    return ( t_cos + t_sin );
}


void        calculate_kochin_pert_coeffs(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            int             freq_pos,
                                            SimulationData* sim_data
                                        )
{
    // Define local variables
    cusfloat    ang_freq    = input->angfreqs[freq_pos];
    int         idx0        = 0;
    int         idx1        = 0;
    int         idx2        = 0;

    // Get current wave number
    cusfloat            wave_num    = w2k( 
                                            ang_freq, 
                                            input->water_depth, 
                                            input->grav_acc 
                                        );
    
    // Define Kochin object interface
    KochinInterface*    kochin      = new KochinInterface(
                                                            mesh_gp->source_nodes[0],
                                                            wave_num,
                                                            input->water_depth,
                                                            0,
                                                            false
                                                        );
    // Calcualte perturbation potential
    cuscomplex*         pert_src    = generate_empty_vector<cuscomplex>( input->heads_np * mesh_gp->source_nodes_tnp );

    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            for ( int k=mesh_gp->source_nodes_cnp[j]; k<mesh_gp->source_nodes_cnp[j+1]; k++ )
            {
                // Update index
                idx0 = i * mesh_gp->source_nodes_tnp + k;

                // Add radiation potential
                for ( int r=0; r<input->dofs_np; r++ )
                {
                    idx1            =  ( 
                                            r * mesh_gp->source_nodes_tnp 
                                            + 
                                            k
                                        );
                    idx2            =  (
                                            i * ( mesh_gp->meshes_np * input->dofs_np )
                                            +
                                            j * input->dofs_np 
                                            + 
                                            r
                                        );
                    
                    pert_src[idx0]  += cuscomplex( 0.0, -1.0 ) * ang_freq * sim_data->raos[idx2] * sim_data->intensities[idx1];
                }

                // Add diffraction potential
                idx1            = ( input->dofs_np + i ) * mesh_gp->source_nodes_tnp + k;
                pert_src[idx0] += sim_data->intensities[idx1];

            }
        }
    }
    
    // Loop over DOFs in order to calculate the kochin raddiation coefficients
    for ( int i=0; i<input->heads_np; i++ )
    {
        idx0    = i * mesh_gp->source_nodes_tnp;
        idx1    =  i * ( mesh_gp->meshes_np * input->kochin_np );

        calculate_kochin_coefficients(
                                        input,
                                        mesh_gp,
                                        kochin,
                                        &(pert_src[idx0]),
                                        &(sim_data->mdrift_kochin_pert_cos[idx1]),
                                        &(sim_data->mdrift_kochin_pert_sin[idx1])
                                    );
    }

    // Delete heap memory
    delete kochin;
    mkl_free( pert_src );

}


void        calculate_kochin_rad_coeffs(
                                            Input*          input,
                                            MeshGroup*      mesh_gp,
                                            int             freq_pos,
                                            SimulationData* sim_data
                                        )
{
    // Define local variables
    int idx0    = 0;
    int idx1    = 0;

    // Get current wave number
    cusfloat    wave_num        = w2k( 
                                            input->angfreqs[freq_pos], 
                                            input->water_depth, 
                                            input->grav_acc 
                                        );
    
    // Define Kochin object interface
    KochinInterface*    kochin  = new KochinInterface(
                                                        mesh_gp->source_nodes[0],
                                                        wave_num,
                                                        input->water_depth,
                                                        0,
                                                        false
                                                    );
    
    // Loop over DOFs in order to calculate the kochin raddiation coefficients
    for ( int i=0; i<input->dofs_np; i++ )
    {
        idx0    = i * mesh_gp->source_nodes_tnp;
        idx1    = i * ( mesh_gp->meshes_np * input->kochin_np );

        calculate_kochin_coefficients(
                                        input,
                                        mesh_gp,
                                        kochin,
                                        &(sim_data->intensities[idx0]),
                                        &(sim_data->mdrift_kochin_rad_cos[idx1]),
                                        &(sim_data->mdrift_kochin_rad_sin[idx1])
                                    );
    }

    // Delete heap memory
    delete kochin;
}