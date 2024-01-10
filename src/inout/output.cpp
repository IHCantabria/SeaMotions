
// Include general usage libraries
#include <filesystem>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "output.hpp"
#include "../tools.hpp"
#include "../version.hpp"


Output::Output( 
                                            Input*      input
                )
{
    // Storage pointer to the input class instance
    this->_input = input;

    // Generate output results file path
    std::filesystem::path _case_fopath( input->case_fopath );
    std::filesystem::path _results_finame( std::string( "results.hydb.h5" ) );
    std::filesystem::path _results_fipath = _case_fopath / _results_finame;

    this->_results_fipath = _results_fipath.string( );


    /******************************************************/
    /************ Define datasets dimensions **************/
    /******************************************************/

    // Storage frequencies and headings datasets dimensions
    this->_ds_fh[0]     = input->angfreqs_np;

    // Storage structural mass and hydrostatic stiffness datasets dimensions
    this->_ds_mh[0]     = input->bodies_np;
    this->_ds_mh[1]     = input->dofs_np;
    this->_ds_mh[2]     = input->dofs_np;

    // Storage hydromechanics coefficients datasets dimensions
    this->_ds_hm[0]    = input->bodies_np;
    this->_ds_hm[1]    = input->bodies_np;
    this->_ds_hm[2]    = input->angfreqs_np;
    this->_ds_hm[3]    = input->dofs_np;
    this->_ds_hm[4]    = input->dofs_np;

    // Storage QTF datasets dimensions for the ith body
    this->_ds_qf[0]   = input->bodies_np;
    this->_ds_qf[1]   = input->heads_np;
    this->_ds_qf[2]   = input->heads_np;
    this->_ds_qf[3]   = input->angfreqs_np;
    this->_ds_qf[4]   = input->angfreqs_np;
    this->_ds_qf[5]   = input->dofs_np;

    // Storage mean drift dataset dimensions for the ith body
    this->_ds_wx[0]    = input->heads_np;
    this->_ds_wx[1]    = input->bodies_np;
    this->_ds_wx[2]    = input->angfreqs_np;
    this->_ds_wx[3]    = input->dofs_np;

    /******************************************************/
    /**** Create datasets for the required output data ****/
    /******************************************************/

    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_TRUNC );

    CREATE_ATTRIBUTE( fid, "version_major", VERSION_MAJOR, H5::PredType::NATIVE_INT );
    CREATE_ATTRIBUTE( fid, "version_minor", VERSION_MINOR, H5::PredType::NATIVE_INT );
    CREATE_ATTRIBUTE( fid, "version_patch", VERSION_PATCH, H5::PredType::NATIVE_INT );

    if ( input->out_diffrac )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_DIFFRAC_MAG,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_DIFFRAC_PHA,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );
    }

    if ( input->out_fk )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_FK_MAG,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_FK_PHA,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );
    }

    if ( input->out_hydmech )
    {
        // Create dataset for added mass
        CREATE_DATASET( 
                            fid,
                            _DN_ADDED_MASS,
                            _DS_HM_NP,
                            this->_ds_hm,
                            cusfloat_h5
                        );

        // Create dataset for radiation damping
        CREATE_DATASET( 
                            fid,
                            _DN_DAMPING_RAD,
                            _DS_HM_NP,
                            this->_ds_hm,
                            cusfloat_h5
                        );

    }

    if ( input->out_hydstiff )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_HYDSTIFF,
                            _DS_MH_NP,
                            this->_ds_mh,
                            cusfloat_h5
                        );
    }

    // Create dataset for the frequencies set
    CREATE_DATASET( 
                        fid,
                        _DN_FREQS,
                        _DS_HF_NP,
                        this->_ds_fh,
                        cusfloat_h5
                    );

    // Create dataset for the headings set
    CREATE_DATASET( 
                        fid,
                        _DN_HEADS,
                        _DS_HF_NP,
                        this->_ds_fh,
                        cusfloat_h5
                    );

    // Create dataset for mean drift forces
    if ( input->is_calc_mdrift )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_MDRIFT_MAG,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );
        
        CREATE_DATASET( 
                            fid,
                            _DN_MDRIFT_PHA,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );

        if ( this->_input->out_qtf_comp )
        {
            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_WL_MAG,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_WL_PHA,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_BERN_MAG,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_BERN_PHA,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_ACC_MAG,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_ACC_PHA,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_MOM_MAG,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_MDRIFT_MOM_PHA,
                                _DS_WX_NP,
                                this->_ds_wx,
                                cusfloat_h5
                            );
        }
    }
    
    // Create dataset for QTF frequency difference
    if ( input->out_qtf )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_QTF_DIFF_MAG,
                            _DS_QF_NP,
                            this->_ds_qf,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_QTF_DIFF_PHA,
                            _DS_QF_NP,
                            this->_ds_qf,
                            cusfloat_h5
                        );

        if ( this->_input->out_qtf_comp )
        {
            // Set ouput for Acceleration component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_ACC_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_ACC_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            // Set output for Bernouilly component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_BERN_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_BERN_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            // Set output for Momentum component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_MOM_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_MOM_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            // Set output for WL component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_WL_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_DIFF_WL_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
        }

        // Create dataset for QTF frequency summation
        CREATE_DATASET( 
                            fid,
                            _DN_QTF_SUM_MAG,
                            _DS_QF_NP,
                            this->_ds_qf,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_QTF_SUM_PHA,
                            _DS_QF_NP,
                            this->_ds_qf,
                            cusfloat_h5
                        );

        if ( this->_input->out_qtf_comp )
        {
            // Set ouput for Acceleration component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_ACC_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_ACC_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            // Set output for Bernouilly component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_BERN_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_BERN_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            // Set output for Momentum component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_MOM_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_MOM_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
            // Set output for WL component
            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_WL_MAG,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );

            CREATE_DATASET( 
                                fid,
                                _DN_QTF_SUM_WL_PHA,
                                _DS_QF_NP,
                                this->_ds_qf,
                                cusfloat_h5
                            );
            
        }

    }

    // Create dataset for RAOs
    if ( input->out_raos )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_RAO_MAG,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_RAO_PHA,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );
        
    }

    // Create dataset for structural mass
    if ( input->out_struct_mass )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_STRUCT_MASS,
                            _DS_MH_NP,
                            this->_ds_mh,
                            cusfloat_h5
                        );
    }

    // Create dataset for wave exciting forces
    if ( input->out_wex )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_WEX_MAG,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_WEX_PHA,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );

    }

    // Close file unit
    fid.close( );

}


void    Output::save_frequencies(
                                            cusfloat*   freqs
                                )
{
    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Storage data
    hsize_t offset_freqs[_DS_HF_NP] = { 0 };
    SAVE_DATASET_CHUNK(
                            fid,
                            _DN_FREQS,
                            _DS_HF_NP,
                            this->_ds_fh,
                            this->_ds_fh,
                            offset_freqs,
                            freqs,
                            cusfloat_h5
                        );

    // Close file unit
    fid.close( );
}


void    Output::save_headings(
                                            cusfloat*   heads
                            )
{
    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Storage data
    hsize_t offset_heads[_DS_HF_NP] = { 0 };
    SAVE_DATASET_CHUNK(
                            fid,
                            _DN_HEADS,
                            _DS_HF_NP,
                            this->_ds_fh,
                            this->_ds_fh,
                            offset_heads,
                            heads,
                            cusfloat_h5
                        );

    // Close file unit
    fid.close( );
}


void    Output::save_hydromechanics_format( 
                                            int         freq_index,
                                            std::string channel_name,
                                            cusfloat*   channel_data
                                            )
{
    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Allocate space for body data chunks
    cusfloat*   amb = generate_empty_vector<cusfloat>( pow2s( this->_input->dofs_np ) );

    // Storage data into disk
    hsize_t _ds_hm_ch[_DS_HM_NP]    = { 1, 1, 1, this->_input->dofs_np, this->_input->dofs_np };
    hsize_t offset[_DS_HM_NP]       = { 0, 0, freq_index, 0, 0 };

    int index_global                = 0;
    int index_local                 = 0;

    for ( int i=0; i<this->_input->bodies_np; i++ )
    {
        for ( int j=0; j<this->_input->bodies_np; j++ )
        {
            // Define body data chunk
            for ( int k=0; k<this->_input->dofs_np; k++ )
            {
                for ( int m=0; m<this->_input->dofs_np; m++ )
                {
                    // Define indexes to access to the required
                    // memory allocation
                    index_global    = (
                                            i * this->_input->bodies_np * pow2s( this->_input->dofs_np )
                                            +
                                            j * this->_input->dofs_np
                                            +
                                            k * this->_input->bodies_np * this->_input->dofs_np
                                            +
                                            m
                                        );
                    index_local     = k * this->_input->dofs_np+m;

                    // Copy data from global matrixes to local body matrixes
                    amb[index_local] = channel_data[index_global];
                }
            }

            // Set dataset chunk offset
            offset[0] = i;
            offset[1] = j;

            // Storage data
            SAVE_DATASET_CHUNK(
                                    fid,
                                    channel_name.c_str( ),
                                    _DS_HM_NP,
                                    this->_ds_hm,
                                    _ds_hm_ch,
                                    offset,
                                    amb,
                                    cusfloat_h5
                                );
        }
    }

    // Close file unit
    fid.close( );

    // Deallocate heap memory
    mkl_free( amb );
}


void    Output::save_hydstiffness(
                                            Hydrostatics** hydrostatics
                                )
{
    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Loop over bodies to storage hydrostatic stiffness data matrix
    hsize_t offset[_DS_MH_NP]       = { 0, 0, 0 };
    hsize_t _ds_mh_ch[_DS_MH_NP]    = { 1, this->_input->dofs_np, this->_input->dofs_np };

    for ( int i=0; i<this->_input->bodies_np; i++ )
    {
        offset[0] = i;
        SAVE_DATASET_CHUNK(
                                fid,
                                _DN_HYDSTIFF,
                                _DS_MH_NP,
                                this->_ds_mh,
                                _ds_mh_ch,
                                offset,
                                hydrostatics[i]->hydstiffmat,
                                cusfloat_h5
                            );
    }

    // Close file unit
    fid.close( );
}


void    Output::save_structural_mass( 
                                            void 
                                    )
{
    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Allocate space for the ith body 
    // structural mass matrix
    cusfloat*   body_mass    = generate_empty_vector<cusfloat>( pow2s( this->_input->dofs_np ) );

    // Storage structural mass matrix for all the bodies 
    // in the simulation
    hsize_t offset[_DS_MH_NP]       = { 0, 0, 0 };
    hsize_t _ds_mh_ch[_DS_MH_NP]    = { 1, this->_input->dofs_np, this->_input->dofs_np };

    for ( int i=0; i<this->_input->bodies_np; i++ )
    {
        // Clear structural mass matrix to delete spurious
        // data from the previous body if any
        clear_vector( pow2s( this->_input->dofs_np ), body_mass );

        // Create ith body mass matrix
        // Define body mass matrix
        body_mass[0]    = this->_input->bodies[i]->mass;  // Surge
        body_mass[7]    = this->_input->bodies[i]->mass;  // Sway
        body_mass[14]   = this->_input->bodies[i]->mass; // Heave

        if ( this->_input->bodies[i]->interia_by_rad )
        {
            body_mass[21] = this->_input->bodies[i]->mass * pow2s( this->_input->bodies[i]->rad_inertia[0] ); // Roll
            body_mass[28] = this->_input->bodies[i]->mass * pow2s( this->_input->bodies[i]->rad_inertia[1] ); // Pitch
            body_mass[35] = this->_input->bodies[i]->mass * pow2s( this->_input->bodies[i]->rad_inertia[2] ); // Yaw
        }
        else
        {
            body_mass[21] = this->_input->bodies[i]->inertia[0]; // Roll
            body_mass[22] = this->_input->bodies[i]->inertia[1]; // Roll - Pitch
            body_mass[23] = this->_input->bodies[i]->inertia[2]; // Roll - Yaw
            body_mass[27] = this->_input->bodies[i]->inertia[1]; // Pitch - Roll
            body_mass[28] = this->_input->bodies[i]->inertia[3]; // Pitch
            body_mass[29] = this->_input->bodies[i]->inertia[4]; // Pitch - Yaw
            body_mass[33] = this->_input->bodies[i]->inertia[2]; // Yaw - Roll
            body_mass[34] = this->_input->bodies[i]->inertia[4]; // Yaw - Pitch-
            body_mass[35] = this->_input->bodies[i]->inertia[5]; // Yaw
        }

        // Storage data
        offset[0] = i;
        SAVE_DATASET_CHUNK(
                                fid,
                                _DN_STRUCT_MASS,
                                _DS_MH_NP,
                                this->_ds_mh,
                                _ds_mh_ch,
                                offset,
                                body_mass,
                                cusfloat_h5
                            );

    }

    // Delete local heap memory
    mkl_free( body_mass );

    // Close file unit
    fid.close( );
}


void    Output::save_qtf_format(
                                            std::string channel_name,
                                            cuscomplex* forces
                            )
{
    // Get local handle to input object
    Input* input    = this->_input;

    // Define datasets name
    std::string _dn_mag = channel_name + std::string( "_mag" );
    std::string _dn_pha = channel_name + std::string( "_pha" );

    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Convert to magnitude and phase format. Also, the data is re-ordered in the format: NB x NH x NH x NF x NF x NDOF
    int         force_np    =  (
                                    input->bodies_np 
                                    *
                                    pow2s( input->heads_np )
                                    *
                                    pow2s( input->angfreqs_np )
                                    *
                                    input->dofs_np
                                );
    cusfloat* data_mag  = generate_empty_vector<cusfloat>( force_np );
    cusfloat* data_pha  = generate_empty_vector<cusfloat>( force_np );

    int id_old  = 0;
    int id_new  = 0;

    for ( int ib=0; ib<this->_input->bodies_np; ib++ )
    {
        for ( int ih1=0; ih1<this->_input->heads_np; ih1++ )
        {
            for ( int ih2=0; ih2<this->_input->heads_np; ih2++ )
            {
                for ( int ifr1=0; ifr1<this->_input->angfreqs_np; ifr1++ )
                {
                    for ( int ifr2=0; ifr2<this->_input->angfreqs_np; ifr2++ )
                    {
                        for ( int id=0; id<this->_input->dofs_np; id++ )
                        {
                            // Define matrix indexes
                            id_old      =  (
                                                ih1  * ( input->heads_np * input->bodies_np * pow2s( input->angfreqs_np ) * input->dofs_np )
                                                +
                                                ih2  * ( input->bodies_np * pow2s( input->angfreqs_np ) * input->dofs_np )
                                                +
                                                ib   * ( pow2s( input->angfreqs_np ) * input->dofs_np )
                                                +
                                                ifr1 * ( input->angfreqs_np * input->dofs_np )
                                                +
                                                ifr2 *  input->dofs_np
                                                +
                                                id
                                            );

                            id_new      =   (
                                                ib   * ( pow2s( input->heads_np ) * pow2s( input->angfreqs_np ) * input->dofs_np )
                                                +
                                                ih1  * ( input->heads_np * pow2s( input->angfreqs_np ) * input->dofs_np )
                                                +
                                                ih2  * ( pow2s( input->angfreqs_np ) * input->dofs_np )
                                                +
                                                ifr1 * ( input->angfreqs_np * input->dofs_np )
                                                +
                                                ifr2 * input->dofs_np
                                                +
                                                id
                                            );

                            // Re-order data
                            data_mag[id_new]    = std::abs( forces[id_old] );
                            data_pha[id_new]    = std::atan2( forces[id_old].imag( ), forces[id_old].real( ) );

                        }
                    }
                }
            }
        }
    }

    // Storage input data into disk
    hsize_t _ds_qf_ch[_DS_QF_NP]    = { input->bodies_np, input->heads_np, input->heads_np, input->angfreqs_np, input->angfreqs_np, input->dofs_np };
    hsize_t offset[_DS_QF_NP]       = { 0, 0, 0, 0, 0, 0 };

    SAVE_DATASET_CHUNK(
                            fid,
                            _dn_mag.c_str( ),
                            _DS_QF_NP,
                            this->_ds_qf,
                            _ds_qf_ch,
                            offset,
                            data_mag,
                            cusfloat_h5
                        );

    SAVE_DATASET_CHUNK(
                            fid,
                            _dn_pha.c_str( ),
                            _DS_QF_NP,
                            this->_ds_qf,
                            _ds_qf_ch,
                            offset,
                            data_pha,
                            cusfloat_h5
                        );

    // Close file unit
    fid.close( );

    // Deallocate heap memory used for the current function
    mkl_free( data_mag );
    mkl_free( data_pha );
}


void    Output::save_wave_exciting_format(
                                            int         freq_index,
                                            std::string channel_name,
                                            cuscomplex* forces
                                        )
{
    // Define datasets name
    std::string _dn_mag = channel_name + std::string( "_mag" );
    std::string _dn_pha = channel_name + std::string( "_pha" );

    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_RDWR );

    // Allocate space for ith body data
    cusfloat* data_mag = generate_empty_vector<cusfloat>( this->_input->dofs_np );
    cusfloat* data_pha = generate_empty_vector<cusfloat>( this->_input->dofs_np );

    // Storage input data into disk
    hsize_t _ds_wx_ch[_DS_WX_NP]    = { 1, 1, 1, this->_input->dofs_np };
    int     index                   = 0;
    hsize_t offset[_DS_WX_NP]       = { 0, 0, freq_index, 0 };

    for ( int i=0; i<this->_input->heads_np; i++ )
    {
        for ( int j=0; j<this->_input->bodies_np; j++ )
        {
            // Process data to have the format of magnitude and phase
            for ( int k=0; k<this->_input->dofs_np; k++ )
            {
                index       = (
                                    i * this->_input->bodies_np * this->_input->dofs_np
                                    +
                                    j * this->_input->dofs_np
                                    +
                                    k
                                );
                data_mag[k] = std::abs( forces[index] );
                data_pha[k] = std::atan2( forces[index].imag( ), forces[index].real( ) );
            }

            // Storage data
            offset[0] = i;
            offset[1] = j;

            SAVE_DATASET_CHUNK(
                                    fid,
                                    _dn_mag.c_str( ),
                                    _DS_WX_NP,
                                    this->_ds_wx,
                                    _ds_wx_ch,
                                    offset,
                                    data_mag,
                                    cusfloat_h5
                                );

            SAVE_DATASET_CHUNK(
                                    fid,
                                    _dn_pha.c_str( ),
                                    _DS_WX_NP,
                                    this->_ds_wx,
                                    _ds_wx_ch,
                                    offset,
                                    data_pha,
                                    cusfloat_h5
                                );

        }
    }

    // Close file unit
    fid.close( );

    // Deallocate heap memory used for the current function
    mkl_free( data_mag );
    mkl_free( data_pha );
}