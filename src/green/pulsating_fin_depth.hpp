
#ifndef __pulsating_hpp
#define __pulsating_hpp

#include "../config.hpp"
#include "../waves.hpp"


cuscomplex john_series(cusfloat R, cusfloat z, cusfloat zeta, cusfloat h,
                        WaveDispersionData &wave_data);

#endif