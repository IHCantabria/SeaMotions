
#ifndef __source_hpp
#define __source_hpp

#include "../config.hpp"
#include "../containers/containers.hpp"


void    calculate_source_monopole_potential(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat&       phi
                                            );


void    calculate_source_monopole_velocity(
                                            PanelGeom*      panel, 
                                            cusfloat        r0, 
                                            cusfloat*       field_point,
                                            cusfloat*       velocity
                                        );


void    calculate_source_potential_newman(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag, 
                                            int             multipole_flag, 
                                            cusfloat&       phi
                                        );


void    calculate_source_velocity_newman(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag,
                                            int             multipole_flag, 
                                            cusfloat        *velocity
                                        );


void    calculate_source_newman(
                                            PanelGeom*      panel, 
                                            cusfloat*       field_point, 
                                            int             fp_local_flag,
                                            int             multipole_flag, 
                                            cusfloat        *velocity,
                                            cusfloat&       phi
                                        );


void    calculate_source_potential_hess(
                                            PanelGeom*      panel,
                                            cusfloat*       field_point,
                                            int             fp_local_flag,
                                            int             multipole_flag,
                                            cusfloat        &phi
                                        );


void    calculate_source_velocity_hess(
                                            PanelGeom*      panel,
                                            cusfloat*       field_point, 
                                            int             fp_local_flag, 
                                            int             multipole_flag,
                                            cusfloat        *velocity
                                        );

#endif