
#ifndef __hmf_interface_hpp
#define __hmf_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../containers/panel_geom.hpp"
#include "grf_interface.hpp"
#include "gwf_interface.hpp"
#include "../inout/input.hpp"


struct HMFInterface
{
private:
    // Define class attributes
    cusfloat        _ang_freq               = 0.0;
    int             _dof_j                  = 0;
    int             _end_index              = 0;
    GRFInterface*   _green_interf_steady    = nullptr;
    GWFInterface*   _green_interf_wave      = nullptr;
    Input*          _input                  = nullptr;
    int             _offset_index           = 0;
    PanelGeom*      _panel                  = nullptr;
    SourceNode**    _source_nodes           = nullptr;
    cuscomplex*     _source_values          = nullptr;
    int             _start_index            = 0;

public:
    // Define class constructor and destructor
    HMFInterface(
                                    SourceNode**    source_nodes,
                                    cuscomplex*     source_values,
                                    PanelGeom*      panel,
                                    int             offset_index,
                                    int             start_index,
                                    int             end_index,
                                    int             dof_j,
                                    cusfloat        ang_freq,
                                    Input*          input_in
                );

    ~HMFInterface(
                                    void
                );

    // Define class methods
    void        set_ang_freq(
                                    cusfloat ang_freq
                            );
    
    cuscomplex  operator()(
                                    cusfloat xi,
                                    cusfloat eta,
                                    cusfloat x,
                                    cusfloat y,
                                    cusfloat z
                            );

    void        set_start_index_i(
                                    int offset_index,
                                    int start_index,
                                    int end_index
                                );
    
    void        set_dof_j(
                                    int dof_j
                        );

    void        set_panel(
                                    PanelGeom* panel
                        );

    void        set_source_values(
                                    cuscomplex* source_values
                                );
};

#endif