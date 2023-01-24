
#ifndef __pulsating_hpp
#define __pulsating_hpp

cuscomplex john_series(cusfloat R, cusfloat z, cusfloat zeta, cusfloat h, cusfloat nu,
                    cusfloat k0, int num_kn, cusfloat* kn);

#endif