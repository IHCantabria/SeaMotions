
#ifndef __gwfdx_interface_hpp
#define __gwfdx_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "gwfdn_interface.hpp"
#include "../waves.hpp"


struct GWFDxInterface: public GWFDnInterface
{
public:
    // Define constructors and destructors
    GWFDxInterface( 
                                SourceNode* source_i,
                                SourceNode* source_j,
                                cusfloat    ang_freq,
                                cusfloat    water_depth,
                                cusfloat    grav_acc
                );

    ~GWFDxInterface(  
                                void 
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