
#ifndef __gwfdn_interface_hpp
#define __gwfdn_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "../waves.hpp"


struct GWFDnInterface
{
protected:
    // Define protected attributes
    cusfloat            _grav_acc           = 0.0;
    cusfloat            _field_point_j[3]   = { 0.0, 0.0, 0.0 };
    SourceNode*         _source_i           = nullptr;
    SourceNode*         _source_j           = nullptr;
    IntegralsDb*        _integrals_db       = nullptr;
    cusfloat            _water_depth        = 0.0;
    WaveDispersionData* _wave_data          = nullptr;

    // Define protected methods
    void    _clear_heap(
                                    void
                        );

    void    _initialize(
                                    cusfloat    ang_freq
                        );

public:
    // Define constructors and destructors
    GWFDnInterface( )   = default;

    GWFDnInterface( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    ang_freq,
                                    cusfloat    water_depth,
                                    cusfloat    grav_acc
                );

    virtual ~GWFDnInterface(  
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

    void        set_ang_freq(
                                    cusfloat    ang_freq
                            );

    void        set_field_point(
                                    cusfloat*   fp
                                );

    void        set_source_i(
                                    SourceNode* source_node
                            );

    void        set_source_j(
                                    SourceNode* source_node
                            );

};

#endif