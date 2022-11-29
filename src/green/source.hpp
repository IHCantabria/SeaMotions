
#ifndef __source_hpp
#define __source_hpp

#include "../config.hpp"
#include "../containers.hpp"


void calculate_source_potential(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat &phi);
void calculate_source_velocity(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag, cusfloat (&velocity)[3]);

#endif