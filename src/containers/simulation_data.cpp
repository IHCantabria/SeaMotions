
// Include local modules
#include "../math/math_tools.hpp"
#include "simulation_data.hpp"


SimulationData::SimulationData(
                                    int         bodies_np,
                                    int         dofs_np,
                                    int         heads_np,
                                    int         rows_local_np,
                                    int         cols_local_np,
                                    int         rows_np,
                                    MpiConfig*  mpi_config_in
                                )
{
    // Storage input arguments into class attributes
    this->hydmech_np        = pow2s( dofs_np * bodies_np );
    this->_mpi_config       = mpi_config_in;  
    this->wave_exc_np       = heads_np * bodies_np * dofs_np;

    // Allocate space variables used in all the processes
    this->added_mass        = generate_empty_vector<cusfloat>( this->hydmech_np );
    this->damping_rad       = generate_empty_vector<cusfloat>( this->hydmech_np );
    this->froude_krylov     = generate_empty_vector<cuscomplex>( this->wave_exc_np );
    this->raos              = generate_empty_vector<cuscomplex>( this->wave_exc_np );
    this->intensities       = generate_empty_vector<cuscomplex>( ( dofs_np + heads_np ) * rows_local_np );
    this->panels_potential  = generate_empty_vector<cuscomplex>( ( dofs_np + heads_np ) * rows_np );
    this->wave_diffrac      = generate_empty_vector<cuscomplex>( this->wave_exc_np );

    // Allocate space for variables used only on root processor
    if ( this->_mpi_config->is_root( ) )
    {
        this->added_mass_p0         = generate_empty_vector<cusfloat>( hydmech_np );
        this->damping_rad_p0        = generate_empty_vector<cusfloat>( hydmech_np );
        this->froude_krylov_p0      = generate_empty_vector<cuscomplex>( wave_exc_np );
        this->hydrostiff_p0         = generate_empty_vector<cusfloat>( hydmech_np );
        this->structural_mass_p0    = generate_empty_vector<cusfloat>( hydmech_np );
        this->wave_diffrac_p0       = generate_empty_vector<cuscomplex>( wave_exc_np );
        this->wave_exc_p0           = generate_empty_vector<cuscomplex>( wave_exc_np );
    }
}


SimulationData::~SimulationData(
                                    void
                                )
{
    mkl_free( this->added_mass );
    mkl_free( this->damping_rad );
    mkl_free( this->froude_krylov );
    mkl_free( this->raos );
    mkl_free( this->intensities );
    mkl_free( this->panels_potential );
    mkl_free( this->sysmat );
    mkl_free( this->sysmat_steady );
    mkl_free( this->wave_diffrac );

    if ( this->_mpi_config->is_root( ) )
    {
        mkl_free( this->added_mass_p0 );
        mkl_free( this->damping_rad_p0 );
        mkl_free( this->froude_krylov_p0 );
        mkl_free( this->hydrostiff_p0 );
        mkl_free( this->structural_mass_p0 );
        mkl_free( this->wave_diffrac_p0 );
        mkl_free( this->wave_exc_p0 );
    }
}