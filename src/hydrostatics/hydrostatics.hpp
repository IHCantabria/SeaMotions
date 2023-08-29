
#ifndef __hydrostatics_hpp
#define __hydrostatics_hpp

// Include local modules
#include "../inout/input.hpp"
#include "../mesh/mesh.hpp"


struct Hydrostatics
{
private:
    // Declare class methods
    void _calculate( 
                        Mesh*       mesh
                    );

public:
    // Declare class attributes
    cusfloat bmx            = 0.0;
    cusfloat bmy            = 0.0;
    cusfloat displacement   = 0.0;
    cusfloat grav_acc       = 0.0;
    cusfloat gmx            = 0.0;
    cusfloat gmy            = 0.0;
    cusfloat kb             = 0.0;
    cusfloat rho_water      = 0.0;
    cusfloat volume         = 0.0;

    // Declare class constructors and destructor
    Hydrostatics( 
                        Mesh*       mesh,
                        cusfloat    rhow,
                        cusfloat    grav_acc
                );

    // Declare class methods
    void print( void );
    
};

#endif