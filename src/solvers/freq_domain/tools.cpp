
// Include local modules
#include "tools.hpp"

#include "../../math/math_interface.hpp"


void    calculate_fields_raddif_lin(
                                        Input*          input,
                                        cuscomplex*     intensities,
                                        MLGCmpx*        field_gp
                                    )
{
    // Loop over RHS to compute all the panel potentials the panels potentials
    cuscomplex  alpha( 1.0, 0.0 );
    cuscomplex  beta( 0.0, 0.0 );
    int         icnx = 1;
    int         icny = 1;
    for ( int i=0; i<( input->dofs_np + input->heads_np ); i++ )
    {
        cblas_gemv<cuscomplex>( 
                                    CblasRowMajor,
                                    CblasNoTrans,
                                    field_gp->sysmat_nrows,
                                    field_gp->sysmat_ncols,
                                    &alpha,
                                    field_gp->sysmat,
                                    field_gp->sysmat_ncols,
                                    &(intensities[i*field_gp->sysmat_nrows+field_gp->start_col]),
                                    icnx,
                                    &beta,
                                    &(field_gp->field_values[i*field_gp->sysmat_nrows]),
                                    icny
                                );
    }
}