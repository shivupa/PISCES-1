//  
//TV:
//This routine is called from larnoldi.cpp
//  this routine performs  y = ift(T*fft( x))+V*x
//  where V is diagonal in postion space and T is diagonal in momentum space for  n-D Hamiltonian and x and y are wavefunctions
// Complex fourier transforms can be replaced by real to real fourier transform.
///* direct FT of x can be done without copying to phi_x*/
//

#include <iostream> 
#include "DVR.h"

   
void DVR::VectorFFT(const double *x, double *y)
{
   m_pkc->c++;
   
  #pragma omp parallel 
   {
     // first apply diagonal potential operator y = V*x
      #pragma omp for
      for (int igr = 0; igr < ngp; igr++)
          y[igr] = v_diag[igr] * x[igr];
  
   //copy the input x vector 
   #pragma omp for
   for (int igr = 0; igr < ngp; igr++)
        phi_x[igr] =  x[igr];
     
    #pragma omp barrier   
  }
     
  //Fourier transform the x vector to momentum space
  fftw_execute(plan_forward);

  //apply T operator on FFT(x) 
  #pragma omp parallel 
  {
     #pragma omp for   
     for(int igr = 0; igr < ngp; igr++) KE_phi_k[igr] = KE_diag[igr] * phi_k[igr];    
     #pragma omp barrier
  }
   
  //backward Fourier transform the T*FFT(x) vector to position space
  fftw_execute(plan_backward);

  #pragma omp parallel 
  {  
  //Normalize the fourier transformation. 
  #pragma omp for
     for(int igr = 0; igr < ngp; igr++ ) KE_phi_x[igr] /= (double)ngp ;
  #pragma omp barrier
     //add y = ift(T*fft( x))+V*x
  #pragma omp for
    for (int igr=0; igr<ngp; igr++) y[igr] +=  KE_phi_x[igr].real() ;      
  } 
 
}
  
