#ifndef VTX_FFTW_H
#define VTX_FFTW_H
#include "defs.h"
#if defined(USE_MKL_DFT)
#include <mkl_dfti.h>
#else
#include "fftw3.h"
#endif
#include <complex>
struct VectorFFT
{

  typedef std::complex<double> Complex;

  Complex *phi_xk;     //Copies vector x to new complex array: x can be directly fourier transformed
  double *phi_x;     //new double
  Complex  *KE_phi_xk; //KE matrix multiplied by FFT(x):T*FFT(x)
  size_t ngp;
  size_t ngp2;
  int n_1dbas[3];

#if defined(USE_MKL_DFT)
  DFTI_DESCRIPTOR_HANDLE  plan_yz;
  DFTI_DESCRIPTOR_HANDLE  plan_x;
#else
  fftw_plan plan_forward;
  fftw_plan plan_backward;
#endif

  explicit VectorFFT(int * n_1dbas);
  ~VectorFFT();
  void apply(const double* restrict x, double* restrict y, const double * restrict v_diag, const double * restrict KE_diag);
};
#endif
