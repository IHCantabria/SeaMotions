
// Include local modules
#include "../../config.hpp"
#include "../../math/math_interface.hpp"
#include "raos.hpp"


void    calculate_raos(
                            Input*          input,
                            cusfloat*       structural_mass,
                            cusfloat*       added_mass,
                            cusfloat*       damping_rad,
                            cusfloat*       hydstiffness,
                            cuscomplex*     wave_diffrac,
                            cuscomplex*     froude_krylov,
                            cusfloat        ang_freq,
                            cuscomplex*     rao
                        )
{
    // Allocate space to hold the system matrix
    cuscomplex* sysmat  = generate_empty_vector<cuscomplex>( pow2s( input->bodies_np * input->dofs_np ) );

    // Clear input rao vector to avoid problems with residual data
    clear_vector( input->bodies_np * input->dofs_np * input->heads_np, rao );

    // Fill in matrix system
    int         index   = 0;
    for ( int i=0; i<input->bodies_np; i++ )
    {
        for ( int j=0; j<input->bodies_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                for ( int m=0; m<input->dofs_np; m++ )
                {
                    // Get structural mass
                    index = (
                                i * input->bodies_np * pow2s( input->dofs_np )
                                +
                                j * input->dofs_np
                                +
                                k * input->bodies_np * input->dofs_np
                                +
                                m
                            );
                    sysmat[index] = cuscomplex( 
                                                    -pow2s( ang_freq ) * ( structural_mass[index] + added_mass[index] ) + hydstiffness[index],
                                                    ang_freq * damping_rad[index]
                                                );
                }
            }
        }
    }

    // Fill in right hand side vector
    for ( int i=0; i<input->heads_np; i++ )
    {
        for ( int j=0; j<input->bodies_np; j++ )
        {
            for ( int k=0; k<input->dofs_np; k++ )
            {
                index = (
                            i * input->bodies_np * input->dofs_np
                            +
                            j * input->dofs_np
                            +
                            k
                        );
                rao[index] = wave_diffrac[index] + froude_krylov[index];
            }
        }
    }

    // Solve system of equations
    int     rows_np = input->bodies_np * input->dofs_np;
    int     info    = 0;
    int*    ipiv    = generate_empty_vector<int>( rows_np );
    gesv<cuscomplex>( 
                        &rows_np,
                        &(input->heads_np),
                        sysmat,
                        &(rows_np),
                        ipiv,
                        rao,
                        &(rows_np),
                        &info
                    );

    mkl_free( ipiv );

    if ( info != 0 )
    {
        std::cerr << "ERROR - Calculating RAO" << std::endl;
        std::cerr << "Error solving system of equations - Info: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    // Deallocate local heap memory
    mkl_free( sysmat );

}