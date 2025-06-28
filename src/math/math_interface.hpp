
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


// Interface for vCos - Vector addition
template<typename T>
inline const auto& lv_cos = vsCos;

template<>
inline const auto& lv_cos<float> = vsCos;

template<>
inline const auto& lv_cos<double> = vdCos;

template<>
inline const auto& lv_cos<std::complex<float>> = vcCos;

template<>
inline const auto& lv_cos<std::complex<double>> = vzCos;


// Interface for vExp - Vector addition
template<typename T>
inline const auto& lv_exp = vsExp;

template<>
inline const auto& lv_exp<float> = vsExp;

template<>
inline const auto& lv_exp<double> = vdExp;

template<>
inline const auto& lv_exp<std::complex<float>> = vcExp;

template<>
inline const auto& lv_exp<std::complex<double>> = vzExp;

// Interface for vLog - Vector addition
template<typename T>
inline const auto& lv_log = vsLn;

template<>
inline const auto& lv_log<float> = vsLn;

template<>
inline const auto& lv_log<double> = vdLn;

template<>
inline const auto& lv_log<std::complex<float>> = vcLn;

template<>
inline const auto& lv_log<std::complex<double>> = vzLn;


// Interface for vSin - Vector addition
template<typename T>
inline const auto& lv_sin = vsSin;

template<>
inline const auto& lv_sin<float> = vsSin;

template<>
inline const auto& lv_sin<double> = vdSin;

template<>
inline const auto& lv_sin<std::complex<float>> = vcSin;

template<>
inline const auto& lv_sin<std::complex<double>> = vzSin;


// Interface for vSub - Vector substraction
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


// Interface for vMult - Vector multiplication
template<typename T>
inline const auto& lv_mult = vsMul;

template<>
inline const auto& lv_mult<float> = vsMul;

template<>
inline const auto& lv_mult<double> = vdMul;

template<>
inline const auto& lv_mult<std::complex<float>> = vcMul;

template<>
inline const auto& lv_mult<std::complex<double>> = vzMul;


// Interface for vDiv - Vector division
template<typename T>
inline const auto& lv_div = vsDiv;

template<>
inline const auto& lv_div<float> = vsDiv;

template<>
inline const auto& lv_div<double> = vdDiv;

template<>
inline const auto& lv_div<std::complex<float>> = vcDiv;

template<>
inline const auto& lv_div<std::complex<double>> = vzDiv;


// Interface for vPow - Vector Power
template<typename T>
inline const auto& lv_pow = vsPow;

template<>
inline const auto& lv_pow<float> = vsPow;

template<>
inline const auto& lv_pow<double> = vdPow;

template<>
inline const auto& lv_pow<std::complex<float>> = vcPow;

template<>
inline const auto& lv_pow<std::complex<double>> = vzPow;


// Interface for vPowx - Scalar Power
template<typename T>
inline const auto& lv_powx = vsPowx;

template<>
inline const auto& lv_powx<float> = vsPowx;

template<>
inline const auto& lv_powx<double> = vdPowx;

template<>
inline const auto& lv_powx<std::complex<float>> = vcPowx;

template<>
inline const auto& lv_powx<std::complex<double>> = vzPowx;


// Interface for vPow3o2 - Computes square root of the cube of each element
template<typename T>
inline const auto& lv_pow3o2 = vsPow3o2;

template<>
inline const auto& lv_pow3o2<float> = vsPow3o2;

template<>
inline const auto& lv_pow3o2<double> = vdPow3o2;


// Interface for Sqrt - Vector Sqrt
template<typename T>
inline const auto& lv_sqrt = vsSqrt;

template<>
inline const auto& lv_sqrt<float> = vsSqrt;

template<>
inline const auto& lv_sqrt<double> = vdSqrt;

template<>
inline const auto& lv_sqrt<std::complex<float>> = vcSqrt;

template<>
inline const auto& lv_sqrt<std::complex<double>> = vzSqrt;


// Interface for Cbrt - Vector Cbrt
template<typename T>
inline const auto& lv_cbrt = vsCbrt;

template<>
inline const auto& lv_cbrt<float> = vsCbrt;

template<>
inline const auto& lv_cbrt<double> = vdCbrt;


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