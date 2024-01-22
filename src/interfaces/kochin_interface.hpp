
#ifndef __kochin_interface_hpp
#define __kochin_interface_hpp

// Include local modules
#include "../config.hpp"
#include "../containers/source_node.hpp"


struct KochinInterface
{
private:
    // Define class private attributes
    cuscomplex  _cf_l           = cuscomplex( 0.0, 0.0 );
    cusfloat    _ep_l           = 0.0;
    bool        _is_cos         = false;
    int         _l_order        = 0;
    SourceNode* _source         = nullptr;
    cusfloat    _water_depth    = 0.0;
    cusfloat    _wave_num       = 0.0;

public:
    // Defincde class constructors and destructor
    KochinInterface(
                        SourceNode* source,
                        cusfloat    wave_number,
                        cusfloat    water_depth,
                        int         l_order,
                        bool        is_cos
                    );

    // Define class public attributes
    cuscomplex  operator()( 
                                cusfloat    xi,
                                cusfloat    eta,
                                cusfloat    x,
                                cusfloat    y,
                                cusfloat    z
                        );
            
    void        set_is_cos(
                                bool        is_cos
                            );

    void        set_l_order( 
                                int         l_order 
                           );

    void        set_source(
                               SourceNode* source
                           );

};

#endif