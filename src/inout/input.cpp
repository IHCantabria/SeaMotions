
// Include general usage libraries
#include <iostream>

// Include local modules
#include "input.hpp"
#include "../math/math_interface.hpp"
#include "../math/math_tools.hpp"


void Input::configure( void )
{
    // Check headings input units
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

    // Check input frequencies units
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

    // Sort frequencies from the lowest to the highest
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

    // Create a vector for the frequencies
    this->freqs = generate_empty_vector<cusfloat>( this->angfreqs_np );
    for ( int i=0; i<this->angfreqs_np; i++ )
    {
        this->freqs[i] = angfreq_to_freq( this->angfreqs[i] );
    }

    // Calculate source nodes
    for ( int i=0; i<this->bodies_np; i++ )
    {
        this->bodies[i]->mesh->define_source_nodes(
                                                        this->poly_order,
                                                        this->bodies[i]->cog
                                                    );
    }

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
    std::cout << "BODY DEFINITION: " << std::endl;
    for ( int i=0; i<this->bodies_np; i++ )
    {
        std::cout << " - Body " << i << ": " << this->bodies_finame[i] << std::endl;
        this->bodies[i]->print( );
    }

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