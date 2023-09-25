
#ifndef __hmf_interface_hpp
#define __hmf_interface_hpp

// Include local modules
#include "../containers/source_node.hpp"
#include "../containers/panel_geom.hpp"
#include "grf_interface.hpp"
#include "gwf_interface.hpp"


struct HMFInterface
{
private:
    // Define class attributes
    cusfloat        _ang_freq               = 0.0;
    cusfloat        _grav_acc               = 0.0;
    int             _dof_j                  = 0;
    GRFInterface*   _green_interf_steady    = nullptr;
    GWFInterface*   _green_interf_wave      = nullptr;
    PanelGeom*      _panel                  = nullptr;
    SourceNode**    _source_nodes           = nullptr;
    int             _source_nodes_np        = 0;
    cuscomplex*     _source_values          = nullptr;
    int             _start_index_i          = 0;
    cusfloat        _water_depth            = 0.0;

public:
    // Define class constructor and destructor
    HMFInterface(
                                    SourceNode**    source_nodes,
                                    cuscomplex*     source_values,
                                    int             source_nodes_np,             
                                    PanelGeom*      panel,
                                    int             start_index_i,
                                    int             dof_j,
                                    cusfloat        ang_freq,
                                    cusfloat        water_depth,
                                    cusfloat        grav_acc
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
                                    int dof_i
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