
#ifndef __gwf_interface_hpp
#define __gwf_interface_hpp

// Include local modules
#include "../containers/panel_geom.hpp"
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "../waves.hpp"


struct GWFInterface
{
private:
    // Define local attributes
    cusfloat            _grav_acc       = 0.0;
    PanelGeom*          _panel_i        = nullptr;
    PanelGeom*          _panel_j        = nullptr;
    IntegralsDb*        _integrals_db   = nullptr;
    cusfloat            _water_depth    = 0.0;
    WaveDispersionData* _wave_data      = nullptr;

public:
    // Define constructors and destructors
    GWFInterface( 
                    PanelGeom*  panel_i,
                    PanelGeom*  panel_j,
                    cusfloat    ang_freq,
                    cusfloat    water_depth,
                    cusfloat    grav_acc
                );

    ~GWFInterface(  
                    void 
                );

    // Define class methods
    cuscomplex  operator()( 
                                cusfloat x,
                                cusfloat y,
                                cusfloat z
                           );

    void        set_ang_freq(
                                cusfloat ang_freq
                            );

    void        set_panel_i(
                                PanelGeom* panel
                            );

    void        set_panel_j(
                                PanelGeom* panel
                            );

};

#endif