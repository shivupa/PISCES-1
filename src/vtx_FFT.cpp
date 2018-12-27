//  
//TV:
//This routine is called from larnoldi.cpp
//  this routine performs  y = ift(T*fft( x))+V*x
//  where V is diagonal in postion space and T is diagonal in momentum space for  n-D Hamiltonian and x and y are wavefunctions
// Complex fourier transforms can be replaced by real to real fourier transform.
///* direct FT of x can be done without copying to phi_x*/
//

#include <iostream> 
#include <mpi.h> 

#include "DVR.h"

#include "timer.hpp"

   
void DVR::VectorFFT_old(const double *x, double *y)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
//  fftw_mpi_init();

   m_pkc->c++;


 //if(rank==0) progress_timer t("VectorFFT_old", verbose);
 progress_timer t("VectorFFT_old", verbose);
   
  #pragma omp parallel 
   {
     // first apply diagonal potential operator y = V*x
      #pragma omp for
      for (int igr = 0; igr < ngp; igr++){
          y[igr] = v_diag[igr] * x[igr];
         // cout<<v_diag[igr]<<endl;
      }
   //copy the input x vector 
   #pragma omp for
   for (int igr = 0; igr < ngp; igr++) {
     //  if(rank==0)cout<<"1 phi_x["<<igr<<"]= "<<phi_x[igr]<<endl;
        phi_x[igr] =  x[igr];
       // phi_x[igr][0] =  x[igr];
       // if(rank==0)cout<<"2 phi_x["<<igr<<"]= "<<x[igr]<<endl;
    }

    int igr = 0;
   // if(rank==0)cout<<" n_1dbas[0] ="<< n_1dbas[0]<<endl;
   // for (int i = 0; i < n_1dbas[0]; i++) 
   //  for (int j = 0; j < n_1dbas[1]; j++) {
   //    for (int k = 0; k < n_1dbas[2]; k++){ 
   //      if(rank==0)cout<< k <<"  "<< j <<"  "<< i <<"  "<<phi_x[igr].real()<<endl;
   //      igr++;
   //    }
   //      if(rank==0)cout<<"  "<<endl;
   //   }

    #pragma omp barrier   
  }
     
  //Fourier transform the x vector to momentum space
  //cout<<"---- 1st fftw forward ---- "<<endl;
  fftw_execute(plan_forward);
 // fftw_mpi_execute_dft(plan_forward, (fftw_complex*)phi_x, (fftw_complex*)phi_k );

  //cout<<"---- 2nd fftw forward ---- "<<endl;
  //if(rank==0)cout<<"ngp = "<<ngp<<endl;
  //apply T operator on FFT(x) 
  #pragma omp parallel 
  {
     #pragma omp for   
     for(int igr = 0; igr < ngp; igr++){
           KE_phi_k[igr] = KE_diag[igr] * phi_k[igr]; 
        //  if(rank==0)cout<<"phi_k["<<igr<<"]= "<<phi_k[igr]<<endl;
          //if(rank==0)cout<<"KE_diag["<<igr<<"]= "<<KE_diag[igr]<<endl;
        //  if(rank==0)cout<<"KE_phi_k["<<igr<<"]= "<<KE_phi_k[igr]<<endl;
     }
     #pragma omp barrier

//    int igr =0;
//    for (int i = 0; i < n_1dbas[0]; i++)
//     for (int j = 0; j < n_1dbas[1]; j++) {
//       for (int k = 0; k < n_1dbas[2]; k++){
//         if(rank==0)cout<< k <<"  "<< j <<"  "<< i <<"  "<<phi_k[igr].real()<<endl;
//         igr++;
//       }
//         if(rank==0)cout<<"  "<<endl;
//      }
  }
   
  //backward Fourier transform the T*FFT(x) vector to position space
  fftw_execute(plan_backward);
 // fftw_mpi_execute_dft(plan_backward, (fftw_complex*)KE_phi_k, (fftw_complex*)KE_phi_x );
  //cout<<"---- fftw backward ---- "<<endl;

  #pragma omp parallel 
  {  
  //Normalize the fourier transformation. 
  #pragma omp for
     for(int igr = 0; igr < ngp; igr++ ) KE_phi_x[igr]/= (double)ngp ;
  #pragma omp barrier
     //add y = ift(T*fft( x))+V*x
  #pragma omp for
    for (int igr=0; igr<ngp; igr++) {
       y[igr] +=  KE_phi_x[igr].real() ;      
      // y[igr] +=  KE_phi_x[igr][0] ;      
     //  if(rank==0)cout<<"y["<<igr<<"]= "<<y[igr]<<endl;
     //  if(rank==0)cout<<"KE_phi_x["<<igr<<"]= "<<KE_phi_x[igr]<<endl;
     //  if(rank==0)cout<<"KE_phi_x["<<igr<<"] real= "<<KE_phi_x[igr].real()<<endl;
     //  if(rank==0)cout<<"KE_phi_x["<<igr<<"] imag= "<<KE_phi_x[igr].imag()<<endl;
    }

  } 
 
}
  
