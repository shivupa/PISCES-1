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
#include "defs.h"

#include "/ihome/crc/install/intel-2017.1.132/intel-mpi-2017.1.132/fftw/3.3.5/include/fftw3-mpi.h"
//#include "fftw3-mpi.h"

//#include "fftw3.h"
#include "timer.hpp"

#include <complex>
using namespace std;
typedef complex<double> Complex;

   
void VectorFFT(const double *x, double *y)
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
//  fftw_mpi_init();

//   m_pkc->c++;

  int verbose=0;
  int ngp = 1728000 ;
 if(rank==0) progress_timer t("VectorFFT", verbose);
   


   double* KE_diag;   ///< final diagonal kinetic energy matrices coordinate x,y,z

  Complex *phi_x;     //Copies vector x to new complex array: x can be directly fourier transformed
  Complex *phi_k;    //forward fourier transform of vector x
  Complex  *KE_phi_k; //KE matrix multiplied by FFT(x):T*FFT(x)
  Complex  *KE_phi_x; //  Backward fourier transform of T*FFT(x)
 
  fftw_plan plan_forward;
  fftw_plan plan_backward;
  int xa1,ya2,za3;



   dVec Vec_v_diag;
   double*  v_diag;  ///< potential energy at the grid points v_diag[product of n_1dbas]
 
    int  n_1dbas[3];
     n_1dbas[0] = 120;
     n_1dbas[1] = 120;
     n_1dbas[2] = 120;
    int nrofpts = n_1dbas[0]*n_1dbas[1]*n_1dbas[2];
    KE_diag = new double[nrofpts];
    phi_x       = new Complex[nrofpts];
    phi_k       = new Complex[nrofpts];
    KE_phi_k    = new Complex[nrofpts];
    KE_phi_x    = new Complex[nrofpts];

    int nthreads = 1;
    nthreads = omp_get_max_threads();
 
 //   printf("FFT Using %u inti threads.\n",fftw_init_threads() );
 //   fftw_plan_with_nthreads(nthreads);
    printf("FFT Using %u OpenMP threads.\n",omp_get_max_threads() );

    //Define the forward and backward FFT plans
         plan_forward   = fftw_plan_dft(3, n_1dbas,(fftw_complex*)phi_x,(fftw_complex*)phi_k, FFTW_FORWARD,  FFTW_ESTIMATE);
         plan_backward  = fftw_plan_dft(3, n_1dbas,(fftw_complex*)KE_phi_k,(fftw_complex*)KE_phi_x, FFTW_BACKWARD, FFTW_ESTIMATE);
    

  #pragma omp parallel 
   {
     // first apply diagonal potential operator y = V*x
      #pragma omp for
      for (int igr = 0; igr < ngp; igr++)
          y[igr] = v_diag[igr] * x[igr];
  
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
         // if(rank==0)cout<<"KE_diag["<<igr<<"]= "<<KE_diag[igr]<<endl;
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
  
