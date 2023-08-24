
// Include general usage libraries
#include <iostream>

// Include local modules
#include "input.hpp"
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
        std::runtime_error( "" );
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
        std::runtime_error( "" );
    }
}


Input::~Input( void )
{
    if ( this->is_bodies )
    {
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