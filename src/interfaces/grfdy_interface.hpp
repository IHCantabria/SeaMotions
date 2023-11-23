
#ifndef __grfdy_interface_hpp
#define __grfdy_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "grfdn_interface.hpp"
#include "../waves.hpp"


struct GRFDyInterface: public GRFDnInterface
{
public:
    // Define constructors and destructors
    GRFDyInterface( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    water_depth
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