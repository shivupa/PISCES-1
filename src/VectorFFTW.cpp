#include <omp.h>
#include <iostream>
#include "timer.hpp"
#include "VectorFFT.hpp"

/** optimizations
 * - reuse FFTW plans and remove allocations per function call
 * - use in-place FFT
 * - remove temporary arrays and assignments
 * - combine y = V*x + ifft(fft(x)*V(g))
 * - reuse phi_xk for KE_phi_xk
 */

VectorFFT::VectorFFT(int *ndim)
{
  int nthreads = omp_get_max_threads();
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);

  n_1dbas[0]=ndim[0];
  n_1dbas[1]=ndim[1];
  n_1dbas[2]=ndim[2];
  int ng2=n_1dbas[0]*n_1dbas[1];
  int ng_h =n_1dbas[2]/2+1;

 std:: cout<<"ng2 and ng_h ="<<ng2<<" "<<ng_h<<"\n";
//  int n_1dbas[3]={ng,ng,ng};
  ngp=static_cast<size_t>(n_1dbas[0])*static_cast<size_t>(n_1dbas[1])*static_cast<size_t>(n_1dbas[2]);
  ngp2=static_cast<size_t>(n_1dbas[0])*static_cast<size_t>(n_1dbas[1])*static_cast<size_t>(n_1dbas[2]/2+1);

//std::cout << "ngp and ngp2 " << ngp<<" "<<ngp2 << "\n";
// unsigned long long int un_ngp=ngp;
// unsigned long long int un_ngp2=ngp2;
  //int nrofpts = n_1dbas[0]*n_1dbas[1]*n_1dbas[2];
  phi_xk    = new Complex[ngp2];
//  phi_xk    = new Complex[ngp];
  phi_x    = new double[ngp]; // new

//  plan_forward   = fftw_plan_dft(3, n_1dbas,(fftw_complex*)phi_xk,(fftw_complex*)phi_xk, FFTW_FORWARD, FFTW_MEASURE); 
//  plan_backward  = fftw_plan_dft(3, n_1dbas,(fftw_complex*)phi_xk,(fftw_complex*)phi_xk, FFTW_BACKWARD, FFTW_MEASURE);

  plan_forward   = fftw_plan_dft_r2c(3, n_1dbas,(double*)phi_x,(fftw_complex*)phi_xk, FFTW_FORWARD); 
  plan_backward  = fftw_plan_dft_c2r(3, n_1dbas,(fftw_complex*)phi_xk,(double*)phi_x, FFTW_BACKWARD);
}


VectorFFT::~VectorFFT()
{
  //cleanup
  fftw_destroy_plan(plan_backward);
  fftw_destroy_plan(plan_forward);
  delete [] phi_xk;
  delete [] phi_x;
}

void VectorFFT::apply(const double *x, double *y, const double *v_diag, const double *KE_diag) 
{
 int verbose=0;
 int ng=n_1dbas[0];
 int ng2=n_1dbas[0]*n_1dbas[1];
 int ng_h =n_1dbas[2]/2+1;

//   printf("start apply ~~~~%d   %d  \n",ng2, ng_h );


 progress_timer t("VectorFFT", verbose);
#pragma omp parallel for simd
  for (size_t igr = 0; igr < ngp; igr++)
  {
//    phi_xk[igr] =  x[igr];
    phi_x[igr] =  x[igr]; //new
//    printf("phi_x [%d]= %f %lf \n",igr, phi_x[igr], x[igr]);
  }

  fftw_execute(plan_forward);

//  printf(" ng2 and ng_h : %d %d\n",ng2, ng_h);
#pragma omp parallel for simd 
   for(size_t i=0; i < ng2; i++) {
     for(size_t j=0; j < ng_h; j++) {
    //  printf("i j : %d %d\n",i,j);
     // printf("phi_xk [%d]= %f \n",ng_h*i+j,phi_xk[ng_h*i+j].real());
      phi_xk[ng_h*i+j] *= KE_diag[ng*i+j];
     //   printf("KE_diag [ %d ] =  %f\n",ng*i+j,KE_diag[ng*i+j] );
      //  printf("phi_xk [ %d : %d]\n",ng_h*i+j, ng*i+j );
     }
   }
/*
#pragma omp parallel for simd 
  for(size_t igr = 0; igr < ngp; igr++){
    //KE_phi_xk[igr] = KE_diag[igr] * phi_xk[igr]; 

    printf("phi_xk [%d]= %lf %lf\n",igr,phi_xk[igr].real(), KE_diag[igr] );
    phi_xk[igr] *= KE_diag[igr];
  }
    printf("============================== \n");
*/

  fftw_execute(plan_backward);

  const double norm=1.0/double(ngp);
//  const double* restrict ke_phi_real=reinterpret_cast<double*>(phi_xk);
  //Normalize the fourier transformation. 
#pragma omp parallel for simd  
  for(size_t igr = 0; igr < ngp; igr++ ) 
  {
    //y[igr]+=norm*ke_phi_real[igr<<1];
    //y[igr] = v_diag[igr] * x[igr]+norm*ke_phi_real[igr<<1];
   // y[igr] = v_diag[igr] * x[igr]+norm*ke_phi_real[2*igr];
    y[igr] = v_diag[igr] * x[igr]+norm*phi_x[igr];
  }
}

