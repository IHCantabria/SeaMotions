
// Include general usage libraries
#include <filesystem>

// Include local modules
#include "output.hpp"
#include "../tools.hpp"


Output::Output( 
                    Input* input
                )
{
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

    // Storage hydromechanics coefficients datasets dimensions
    this->_ds_hm[0]    = input->bodies_np;
    this->_ds_hm[1]    = input->bodies_np;
    this->_ds_hm[2]    = input->angfreqs_np;
    this->_ds_hm[3]    = input->dofs_np;
    this->_ds_hm[4]    = input->dofs_np;

    // Storage QTF datasets dimensions for the ith body
    this->_ds_qf[0]   = input->bodies_np;
    this->_ds_qf[1]   = input->heads_np;
    this->_ds_qf[2]   = input->angfreqs_np;
    this->_ds_qf[3]   = input->angfreqs_np;
    this->_ds_qf[4]   = 2;
    this->_ds_qf[5]   = input->dofs_np;

    // Storage mean drift dataset dimensions for the ith body
    this->_ds_wx[0]    = input->bodies_np;
    this->_ds_wx[1]    = input->heads_np;
    this->_ds_wx[2]    = input->angfreqs_np;
    this->_ds_wx[3]    = input->dofs_np;

    /******************************************************/
    /**** Create datasets for the required output data ****/
    /******************************************************/

    // Open file unit
    H5::H5File fid( this->_results_fipath.c_str( ), H5F_ACC_TRUNC );

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

    // Create dataset for the frequencies set
    CREATE_DATASET( 
                        fid,
                        _DN_FREQS,
                        _DS_HF_NP,
                        this->_ds_fh,
                        cusfloat_h5
                    );

    hsize_t offset_freqs[_DS_HF_NP] = { 0 };
    SAVE_DATASET_CHUNK(
                            fid,
                            _DN_FREQS,
                            _DS_HF_NP,
                            this->_ds_fh,
                            this->_ds_fh,
                            offset_freqs,
                            input->freqs,
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

    hsize_t offset_heads[_DS_HF_NP] = { 0 };
    SAVE_DATASET_CHUNK(
                            fid,
                            _DN_HEADS,
                            _DS_HF_NP,
                            this->_ds_fh,
                            this->_ds_fh,
                            offset_heads,
                            input->heads.data( ),
                            cusfloat_h5
                        );
    
    // Create dataset for mean drift forces
    if ( input->out_mdrift )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_MDRIFT,
                            _DS_WX_NP,
                            this->_ds_wx,
                            cusfloat_h5
                        );
    }
    
    // Create dataset for QTF frequency difference
    if ( input->out_qtf )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_QTF_DIFF,
                            _DS_QF_NP,
                            this->_ds_qf,
                            cusfloat_h5
                        );

        // Create dataset for QTF frequency summation
        CREATE_DATASET( 
                            fid,
                            _DN_QTF_SUM,
                            _DS_QF_NP,
                            this->_ds_qf,
                            cusfloat_h5
                        );

    }

    // Create dataset for RAOs
    if ( input->out_raos )
    {
        CREATE_DATASET( 
                            fid,
                            _DN_RAO_MAG,
                            _DS_HM_NP,
                            this->_ds_hm,
                            cusfloat_h5
                        );

        CREATE_DATASET( 
                            fid,
                            _DN_RAO_PHA,
                            _DS_HM_NP,
                            this->_ds_hm,
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