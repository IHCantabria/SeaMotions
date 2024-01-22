
// Include local modules
#include "kochin.hpp"

#include "../math/integration.hpp"


void    calculate_kochin_coefficients(
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