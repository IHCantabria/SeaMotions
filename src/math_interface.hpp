
#ifndef __math_interface_hpp
#define __math_interface_hpp

// Include general usage libraries
#include <complex>

// Include general usage scientific libraries
#include "mkl.h"


// Interface for vAdd - Vector addition
template<typename T>
inline const auto& lv_add = vsAdd;

template<>
inline const auto& lv_add<float> = vsAdd;

template<>
inline const auto& lv_add<double> = vdAdd;

template<>
inline const auto& lv_add<std::complex<float>> = vcAdd;

template<>
inline const auto& lv_add<std::complex<double>> = vzAdd;


// Interface for cblas_nrm2 - Vector Norm
template<typename T>
inline const auto& cblas_nrm2 = cblas_snrm2;

template<>
inline const auto& cblas_nrm2<float> = cblas_snrm2;

template<>
inline const auto& cblas_nrm2<double> = cblas_dnrm2;

template<>
inline const auto& cblas_nrm2<std::complex<float>> = cblas_scnrm2;

template<>
inline const auto& cblas_nrm2<std::complex<double>> = cblas_dznrm2;


#endif