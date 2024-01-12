
#ifndef __wave_dispersion_base_hpp
#define __wave_dispersion_base_hpp

// Include local modules
#include "../config.hpp"


cusfloat    k2w( 
                                            cusfloat k, 
                                            cusfloat h, 
                                            cusfloat g 
                );

cusfloat    w2k(
                                            cusfloat w, 
                                            cusfloat h, 
                                            cusfloat g
                );

void        w2ki(
                                            cusfloat w, 
                                            cusfloat h, 
                                            cusfloat g, 
                                            int n, 
                                            cusfloat* kn
                );

#endif