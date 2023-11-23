
#ifndef __grfdn_interface_hpp
#define __grfdn_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "../waves.hpp"


struct GRFDnInterface
{
protected:
    // Define local attributes
    cusfloat            _field_point_j[3]   = { 0.0, 0.0, 0.0 };
    SourceNode*         _source_i           = nullptr;
    SourceNode*         _source_j           = nullptr;
    cusfloat            _water_depth        = 0.0;

public:
    // Define constructors and destructors
    GRFDnInterface( )   = default;

    GRFDnInterface( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    water_depth
                );

    virtual ~GRFDnInterface(  
                                    void 
                );

    // Define class methods
    cuscomplex  operator()( 
                                    cusfloat    xi,
                                    cusfloat    eta,
                                    cusfloat    x,
                                    cusfloat    y,
                                    cusfloat    z
                           );

    void        set_field_point_j(
                                    cusfloat* fp
                                );

    void        set_source_i(
                                    SourceNode* source_node
                            );

    void        set_source_j(
                                    SourceNode* source_node
                            );

};

#endif