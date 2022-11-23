
#ifndef __math_interface_hpp
#define __math_interface_hpp

#include <complex>
#include "mkl.h"

// Interface for cblas_nrm2 - Vector Norm
template<typename T>
const auto& cblas_nrm2 = cblas_snrm2;

template<>
const auto& cblas_nrm2<float> = cblas_snrm2;

template<>
const auto& cblas_nrm2<double> = cblas_dnrm2;

template<>
const auto& cblas_nrm2<std::complex<float>> = cblas_scnrm2;

template<>
const auto& cblas_nrm2<std::complex<double>> = cblas_dznrm2;


#endif