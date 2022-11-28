
#ifndef __dipole_hpp
#define __dipole_hpp

#include "../config.hpp"
#include "../containers.hpp"

cusfloat dipole_potential(PanelGeom &panel, cusfloat (&field_point)[3], int fp_local_flag);

#endif