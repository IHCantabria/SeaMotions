
// Include local modules
#include "gwf_interface.hpp"


GWFInterface::GWFInterface( 
                                PanelGeom*  panel_i,
                                PanelGeom*  panel_j,
                                cusfloat    ang_freq,
                                cusfloat    water_depth,
                                cusfloat    grav_acc
                            )
{
    // Storage panels memory address, water depth and 
    // gravitational acceleration to use along the class
    this->_grav_acc     = grav_acc;
    this->_panel_i      = panel_i;
    this->_panel_j      = panel_j;
    this->_water_depth  = water_depth;

    // Calculate wave numbers
    this->_wave_data    = new WaveDispersionData( 
                                                    ang_freq,
                                                    30,
                                                    water_depth,
                                                    grav_acc
                                                );

    // Load integrals database
    this->_integrals_db = new IntegralsDb( );

}


GWFInterface::~GWFInterface( void )
{
    std::cerr << "Calling GWFInterface destructor..." << std::endl;
    delete this->_integrals_db;
    delete this->_wave_data;
}


cuscomplex  GWFInterface::operator()( 
                                        cusfloat x,
                                        cusfloat y,
                                        cusfloat z
                                    )
{
    // Calculate horizontal radius
    cusfloat R = std::sqrt(
                                pow2s( this->_panel_j->center[0] - x )
                                +
                                pow2s( this->_panel_j->center[1] - y )
                            );

    // Calculate Green function derivatives
    cuscomplex  dG_dR   = G_integral_dr(
                                            R,
                                            this->_panel_j->center[2],
                                            z,
                                            this->_water_depth,
                                            *(this->_wave_data),
                                            *(this->_integrals_db)
                                        );
    cuscomplex  dG_dZ   = G_integral_dz(
                                            R,
                                            this->_panel_j->center[2],
                                            z,
                                            this->_water_depth,
                                            *(this->_wave_data),
                                            *(this->_integrals_db)
                                        );

    // Calculate X and Y cartesian coordinates derivatives
    cusfloat    dX      = this->_panel_j->center[0] - x;
    cuscomplex  dG_dX   = dG_dR * dX;
    cusfloat    dY      = this->_panel_j->center[1] - y;
    cuscomplex  dG_dY   = dG_dR * dY;

    // Calculate normal derivate
    cuscomplex  dG_dn   =   (
                                dG_dR * dG_dX * this->_panel_i->normal_vec[0]
                                +
                                dG_dR * dG_dY * this->_panel_i->normal_vec[1]
                                +
                                dG_dZ * this->_panel_i->normal_vec[2]
                            );

    return dG_dn;
}


void    GWFInterface::set_ang_freq(
                                        cusfloat ang_freq
                                    )
{
    // Delete previous wave data instance
    delete this->_wave_data;

    // Create new wave data for the new angular frequency
    this->_wave_data = new WaveDispersionData( 
                                                    ang_freq,
                                                    30,
                                                    this->_water_depth,
                                                    this->_grav_acc
                                                );
}


void    GWFInterface::set_panel_i(
                                    PanelGeom* panel
                                )
{
    this->_panel_i = panel;
}


void    GWFInterface::set_panel_j(
                                    PanelGeom* panel
                                )
{
    this->_panel_j = panel;
}