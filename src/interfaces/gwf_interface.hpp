
#ifndef __gwf_interface_hpp
#define __gwf_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/integrals_db.hpp"
#include "../waves/wave_dispersion_fo.hpp"


struct GWFInterface
{
private:
    // Define class attributes
    cusfloat            _field_point[3] = { 0.0, 0.0, 0.0 };
    cusfloat            _grav_acc       = 0.0;
    IntegralsDb*        _integrals_db   = nullptr;
    SourceNode*         _source         = nullptr;
    cuscomplex          _source_value   = 0.0;
    cusfloat            _water_depth    = 0.0;
    WaveDispersionFO* _wave_data      = nullptr;

public:

    // Define class constructors and destructor
    GWFInterface(
                                SourceNode*         source_in,
                                cuscomplex          source_value_in,
                                cusfloat*           field_point_in,
                                cusfloat            ang_freq,
                                cusfloat            water_depth_in,
                                cusfloat            grav_acc_in
                );

    ~GWFInterface(
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
                                    cusfloat*   field_point
                                );

    void        set_source_i(
                                    SourceNode* source,
                                    cuscomplex  source_val
                           );
};

#endif