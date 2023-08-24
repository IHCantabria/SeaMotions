
// Include general usage libraries
#include <iostream>

// Include local modules
#include "input.hpp"
#include "../math/math_tools.hpp"


void Input::configure( void )
{
    // Check headings input units
    if ( this->heads_units.compare( "deg" ) )
    {
        for ( int i=0; i<this->heads_np; i++ )
        {
            this->heads[i] = deg_to_rad( this->heads[i] );
        }
    }

    // Check input frequencies units
    if ( this->freqs_unit.compare( "period" ) )
    {
        for ( int i=0; i<this->angfreqs_np; i++ )
        {
            this->angfreqs[i] = period_to_angfreq( this->angfreqs[i] );
        }
    }
    else if ( this->freqs_unit.compare( "freqs" ) )
    {
        for ( int i=0; i<this->angfreqs_np; i++ )
        {
            this->angfreqs[i] = freq_to_angfreq( this->angfreqs[i] );
        }
    }
}


void Input::print( void )
{
    std::cout << std::endl;
    std::cout << "BODY DEFINITION: " << std::endl;
    for ( int i=0; i<this->bodies_np; i++ )
    {
        std::cout << " - Body " << i << ": " << this->bodies_finame[i] << std::endl;
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