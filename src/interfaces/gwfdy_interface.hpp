
#ifndef __gwfdy_interface_hpp
#define __gwfdy_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "gwfdn_interface.hpp"
#include "../waves/wave_dispersion_fo.hpp"


struct GWFDyInterface: public GWFDnInterface
{
public:
    // Define constructors and destructors
    GWFDyInterface( 
                                SourceNode* source_i,
                                SourceNode* source_j,
                                cusfloat    ang_freq,
                                cusfloat    water_depth,
                                cusfloat    grav_acc
                );
    
    // Define class methods
    virtual cuscomplex  operator()( 
                                        cusfloat    xi,
                                        cusfloat    eta,
                                        cusfloat    x,
                                        cusfloat    y,
                                        cusfloat    z
                                );

};

#endif