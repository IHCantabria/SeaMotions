
// Include local modules
#include "../../math/math_tools.hpp"
#include "qtf.hpp"


void    calculate_mean_drift(
                                Input*          input,
                                MeshGroup*      mesh_gp,
                                cuscomplex*     rel_we,
                                cusfloat        ang_freq,
                                cuscomplex*     mean_drift
                            )
{
    // Define aux variables to be used along the function
    int         idx0    = 0;
    int         idx1    = 0;
    cuscomplex  int_mod = 0.0;
    PanelGeom*  panel_k = nullptr;

    // Calculate second order force due to relative wave
    // elevation at the WL
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;
            for ( int k=mesh_gp->panels_wl_cnp[j]; k<mesh_gp->panels_wl_cnp[j+1]; j++ )
            {
                panel_k             = mesh_gp->panels_wl[k];
                idx1                = mesh_gp->panels_wl_tnp + k;
                int_mod             = pow2s( rel_we[idx1] ) * panel_k->len_wl;

                for ( int r=0; r<input->dofs_np; r++ )
                {
                    mean_drift[idx0+r] += int_mod * panel_k->normal_vec[r];
                }
            }

            for ( int r=0; r<input->dofs_np; r++ )
            {
                mean_drift[idx0+r] *= - 0.5 * input->grav_acc * input->water_density;
            }
        }
    }

    // Calculate second order force due to the bernouilly contribution
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<mesh_gp->meshes_np; j++ )
        {
            idx0 = i * ( input->dofs_np * input->bodies_np ) + j * input->bodies_np;
            for ( int k=mesh_gp->panels_cnp[j]; k<mesh_gp->panels_cnp[j+1]; j++ )
            {
                // Calculate
            }
        }
    }

}


void    _calculate_velocity_bc(
                                    PanelGeom*  panel,
                                    cuscomplex* vel_bc
                                )
{
    // Calculate velocity due to the incident wave potential
    cuscomplex vel_w( 0.0, 0.0 );

    // Calculate velocity due to hte
    cuscomplex vel_d( 0.0, 0.0 );
    cuscomplex vel_r( 0.0, 0.0 );
}