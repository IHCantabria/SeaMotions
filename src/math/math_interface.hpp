
#ifndef __math_interface_hpp
#define __math_interface_hpp

// Include general usage libraries
#include <complex>

// Include general usage scientific libraries
#include "mkl.h"
#include "mkl_scalapack.h"


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


// Interface for vAdd - Vector addition
template<typename T>
inline const auto& lv_sub = vsSub;

template<>
inline const auto& lv_sub<float> = vsSub;

template<>
inline const auto& lv_sub<double> = vdSub;

template<>
inline const auto& lv_sub<std::complex<float>> = vcSub;

template<>
inline const auto& lv_sub<std::complex<double>> = vzSub;


// Interface for cblas_dot - Dot product
template<typename T>
inline const auto& cblas_dot = cblas_sdot;

template<>
inline const auto& cblas_dot<float> = cblas_sdot;

template<>
inline const auto& cblas_dot<double> = cblas_ddot;

// Interface for cblas_gemm - Matrix-Matrix product
template<typename T>
inline const auto& cblas_gemm = cblas_sgemm;

template<>
inline const auto& cblas_gemm<float> = cblas_sgemm;

template<>
inline const auto& cblas_gemm<double> = cblas_dgemm;

template<>
inline const auto& cblas_gemm<std::complex<float>> = cblas_cgemm;

template<>
inline const auto& cblas_gemm<std::complex<double>> = cblas_zgemm;

// Interface for cblas_gemv - Matrix-vector product
template<typename T>
inline const auto& cblas_gemv = cblas_sgemv;

template<>
inline const auto& cblas_gemv<float> = cblas_sgemv;

template<>
inline const auto& cblas_gemv<double> = cblas_dgemv;

template<>
inline const auto& cblas_gemv<std::complex<float>> = cblas_cgemv;

template<>
inline const auto& cblas_gemv<std::complex<double>> = cblas_zgemv;


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

// Interface for lapack routines
template<typename T>
inline const auto& gesv = sgesv;

template<>
inline const auto& gesv<float> = sgesv;

template<>
inline const auto& gesv<double> = dgesv;

template<>
inline const auto& gesv<std::complex<float>> = cgesv;

template<>
inline const auto& gesv<std::complex<double>> = zgesv;

// Interface for number sorting functions
template<typename T>
inline const auto& lasrt2           = slasrt2;

template<>
inline const auto& lasrt2<float>    = slasrt2;

template<>
inline const auto& lasrt2<double>   = dlasrt2;

// Interface for ScaLapack routines
template<typename T>
inline const auto& pgesv = psgesv;

template<>
inline const auto& pgesv<float> = psgesv;

template<>
inline const auto& pgesv<double> = pdgesv;

template<>
inline const auto& pgesv<std::complex<float>> = pcgesv;

template<>
inline const auto& pgesv<std::complex<double>> = pzgesv;


#endif