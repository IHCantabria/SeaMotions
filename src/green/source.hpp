
#ifndef __source_hpp
#define __source_hpp

#include "../config.hpp"
#include "../containers/containers.hpp"


void calculate_source_potential_newman(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, int multipole_flag, cusfloat &phi);
void calculate_source_velocity_newman(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, int multipole_flag, cusfloat (&velocity)[3]);
void calculate_source_potential_hess(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, int multipole_flag, cusfloat &phi);
void calculate_source_velocity_hess(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, int multipole_flag, cusfloat (&velocity)[3]);

#endif