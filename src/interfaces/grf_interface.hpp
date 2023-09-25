
#ifndef __grf_interface_hpp
#define __grf_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../waves.hpp"


struct GRFInterface
{
private:
    // Define class attributes
    cusfloat*           _field_point    = nullptr;
    SourceNode*         _source         = nullptr;
    cuscomplex          _source_value   = 0.0;
    cusfloat            _water_depth    = 0.0;

public:

    // Define class constructors and destructor
    GRFInterface(
                                SourceNode*         source_in,
                                cuscomplex          source_value_in,
                                cusfloat*           field_point_in,
                                cusfloat            water_depth_in
                );

    ~GRFInterface(
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

    void        set_field_point(
                                    cusfloat*   field_point
                                );

    void        set_source(
                                    SourceNode* source,
                                    cuscomplex  source_val
                           );
};

#endif