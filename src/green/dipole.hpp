
#ifndef __dipole_hpp
#define __dipole_hpp

#include "../config.hpp"
#include "../containers/containers.hpp"

void    calculate_dipole_potential(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag, 
                                            cusfloat&       phi
                                    );


void    calculate_dipole_potential_kernel(
                                            PanelGeom*      panel, 
                                            cusfloat*       node_fieldp_mod,
                                            cusfloat*       node_fieldp_dx,
                                            cusfloat*       node_fieldp_dy, 
                                            cusfloat*       node_fieldp_dz, 
                                            cusfloat*       delta_xi, 
                                            cusfloat*       delta_eta, 
                                            cusfloat&       phi
                                        );


void    calculate_dipole_velocity(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point,
                                            int             fp_local_flag, 
                                            cusfloat        *velocity
                                );


void    calculate_dipole_velocity_kernel(
                                            PanelGeom*      panel, 
                                            cusfloat*       node_fieldp_mod,
                                            cusfloat*       node_fieldp_dx, 
                                            cusfloat*       node_fieldp_dy, 
                                            cusfloat*       node_fieldp_dz,
                                            cusfloat*       delta_xi,
                                            cusfloat*       delta_eta, 
                                            cusfloat*       velocity
                                        );

#endif