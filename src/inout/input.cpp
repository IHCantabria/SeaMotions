
// Include general usage libraries
#include <iostream>

// Include local modules
#include "input.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"


void Input::configure( void )
{
    /**********************************************************/
    /****************** Check input flags *********************/
    /**********************************************************/

    // Check if fast mode is configured properly
    if ( 
            this->poly_order > 0 
            &&
            this->is_fast_solver
        )
    {
        std::cerr << "Using Fast Mode with high order polynomials ";
        std::cerr << "( > 0) is not allowed!" << std::endl;
        throw std::runtime_error( "" );
    }

    // Check if it is necessary to calculate Mean Drift values
    // This option will allow the frequency calculation of the 
    // QTF depedent signals
    this->is_calc_mdrift = ( this->out_mdrift || this->out_qtf );

    // Check if it is necessary to load the free surface 
    // QTF mesh
    this->is_fs_qtf = this->out_qtf_so_model > 0 ? true: false;

    /**********************************************************/
    /************** Check headings input units ****************/
    /**********************************************************/
    if ( this->heads_units.compare( "deg" ) == 0 )
    {
        for ( int i=0; i<this->heads_np; i++ )
        {
            this->heads[i] = deg_to_rad( this->heads[i] );
        }
    }
    else if ( this->heads_units.compare( "rad" ) != 0 )
    {
        std::cout << std::endl;
        std::cout << "ERROR - INPUT:" << std::endl;
        std::cout << "HeadUnits: " << this->heads_units << " is not a valid parameter." << std::endl;
        std::cout << "Valid heading units are: deg | rad." << std::endl;
        std::cout << std::endl;
        throw std::runtime_error( "" );
    }

    /**********************************************************/
    /************** Check input frequencies units *************/
    /**********************************************************/
    if ( this->freqs_unit.compare( "period" ) == 0 )
    {
        for ( int i=0; i<this->angfreqs_np; i++ )
        {
            this->angfreqs[i] = period_to_angfreq( this->angfreqs[i] );
        }
    }
    else if ( this->freqs_unit.compare( "freq" ) == 0 )
    {
        for ( int i=0; i<this->angfreqs_np; i++ )
        {
            this->angfreqs[i] = freq_to_angfreq( this->angfreqs[i] );
        }
    }
    else if ( this->freqs_unit.compare( "angfreq" ) != 0 )
    {
        std::cout << std::endl;
        std::cout << "ERROR - INPUT:" << std::endl;
        std::cout << "FreqUnit: " << this->freqs_unit << " is not a valid parameter." << std::endl;
        std::cout << "Valid frequency units are: period | freq | angfreq." << std::endl;
        std::cout << std::endl;
        throw std::runtime_error( "" );
    }

    /**********************************************************/
    /**** Sort frequencies from the lowest to the highest *****/
    /**********************************************************/
    int* sort_keys  = generate_empty_vector<int>( this->angfreqs_np );
    int  info       = 0;
    lasrt2<cusfloat>( 
                        "I",
                        &this->angfreqs_np,
                        this->angfreqs.data(),
                        sort_keys,
                        &info
                    );

    if ( info != 0 )
    {
        std::cerr << "ERROR - lasrt2" << std::endl;
        std::cerr << "Sort algorithm could not sort the input angular frequencies.";
        std::cerr << " - Error Code: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    mkl_free( sort_keys );

    /**********************************************************/
    /********** Create a vector for the frequencies ***********/
    /**********************************************************/
    this->freqs = generate_empty_vector<cusfloat>( this->angfreqs_np );
    for ( int i=0; i<this->angfreqs_np; i++ )
    {
        this->freqs[i] = angfreq_to_freq( this->angfreqs[i] );
    }

    /**********************************************************/
    /********** Detect points over the free surface ***********/
    /**********************************************************/
    this->is_wl_points = this->is_calc_mdrift;
    if ( 
            this->is_wl_points
        )
    {
        // Detect points over the free surface
        for ( int i=0; i<this->bodies_np; i++ )
        {
            this->bodies[i]->mesh->detect_wl_points( this->wl_det_prec );
        }
    }

    /**********************************************************/
    /**************** Calculate source nodes ******************/
    /**********************************************************/
    for ( int i=0; i<this->bodies_np; i++ )
    {
        this->bodies[i]->mesh->define_source_nodes(
                                                        this->poly_order,
                                                        this->bodies[i]->cog
                                                    );
    }

}


int  Input::gauss_np_factor_1d( void )
{
    int gf = this->gauss_order;
    // if ( !this->is_block_adaption )
    // {
    //     gf = 1;
    // }

    return gf;
}


int  Input::gauss_np_factor_2d( void )
{
    int gf = pow2s( this->gauss_order );
    // if ( !this->is_block_adaption )
    // {
    //     gf = 1;
    // }

    return gf;
}

Input::~Input( void )
{
    if ( this->is_bodies )
    {
        // Delete frequencies containers
        mkl_free( this->freqs );

        // Delete BodyDef object instances
        for ( int i=0; i<this->bodies_np; i++ )
        {
            delete this->bodies[i];
        }
        
        // Delete vector of BodyDef pointers
        delete [] this->bodies;

    }
}


void Input::print( void )
{
    std::cout << std::endl;
    std::cout << "SOLVER CONTROLS: " << std::endl;
    std::cout << " - Block Adaption: " << this->is_block_adaption << std::endl;
    std::cout << " - Fast Solver: " << this->is_fast_solver << std::endl;
    std::cout << " - Gauss Order: " << this->gauss_order << std::endl;
    std::cout << " - GFDnAbsErr: " << this->gfdn_abs_err << std::endl;
    std::cout << " - GFDnRelErr: " << this->gfdn_rel_err << std::endl;
    std::cout << " - KochinNC: " << this->kochin_np << std::endl;
    std::cout << " - LogSingAna: " << this->is_log_sin_ana << std::endl;
    std::cout << " - PolyOrder: " << this->poly_order << std::endl;
    std::cout << " - PotAbsErr: " << this->pot_abs_err << std::endl;
    std::cout << " - PotRelErr: " << this->pot_rel_err << std::endl;
    std::cout << " - PressAbsErr: " << this->press_abs_err << std::endl;
    std::cout << " - PressRelErr: " << this->press_rel_err << std::endl;
    std::cout << " - QTFSOModel: " << this->out_qtf_so_model << std::endl;
    std::cout << " - WLDetPrec: " << this->wl_det_prec << std::endl;

    std::cout << std::endl;
    std::cout << "BODY DEFINITION: " << std::endl;
    for ( int i=0; i<this->bodies_np; i++ )
    {
        std::cout << " - Body " << i << ": " << this->bodies_finame[i] << std::endl;
        this->bodies[i]->print( );
    }

    std::cout << std::endl;
    std::cout << "OUTPUT CHANNELS: " << std::endl;
    std::cout << " - OutDiffrac: " << this->out_diffrac << std::endl;
    std::cout << " - OutFK: " << this->out_fk << std::endl;
    std::cout << " - OutHydMech: " << this->out_hydmech << std::endl;
    std::cout << " - OutHydStiff: " << this->out_hydstiff << std::endl;
    std::cout << " - OutPress: " << this->out_pressure << std::endl;
    std::cout << " - OutMDrift: " << this->out_mdrift << std::endl;
    std::cout << " - OutMesh: " << this->out_mesh << std::endl;
    std::cout << " - OutQTF: " << this->out_qtf << std::endl;
    std::cout << " - OutQTFComp: " << this->out_qtf_comp << std::endl;
    std::cout << " - OutRAOs: " << this->out_raos << std::endl;
    std::cout << " - OutSources: " << this->out_sources << std::endl;
    std::cout << " - OutStMass: " << this->out_struct_mass << std::endl;
    std::cout << " - OutWex: " << this->out_wex << std::endl;

    std::cout << std::endl;
    std::cout << "SITE CONDITIONS: " << std::endl;
    std::cout << " - Water Density: " << this->water_density << std::endl;
    std::cout << " - Grav. Acc: " << this->grav_acc << std::endl;
    std::cout << " - Water Depth: " << this->water_depth << std::endl;

    std::cout << std::endl;
    std::cout << "HEADINGS: " << std::endl;
    print_vector( this->heads_np, this->heads.data( ), 1, 6 );

    std::cout << std::endl;
    std::cout << "ANGULAR FREQUENCIES: " << std::endl;
    print_vector( this->angfreqs_np, this->angfreqs.data( ), 1, 6 );
    std::cout << std::endl;
}