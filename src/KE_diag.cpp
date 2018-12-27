//TV: Most of the part copied from original sine_dvr.cpp file
// Only Diagonal elements of Kinetic energy operator are calculated
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;

#define SQR(a) ((a)*(a))
#ifndef PI
#define PI 3.14159265358979
#endif


//////////////////////////////
//
//  for FFT
//
//  input:
//    ngp     : no of grid points
//    mu      : mass of the particle
//    x0      : left wall  :  leftmost gridpoint is x1
//    xN+1    : right wall :  rightmost gridpoint is xN
//    ldt     : leading dimension of ekin (if ldt = 0, ekin is stored as a lower triangle)
//    verbose : output level 
//
//  output:
//    x       : grid points (x[0] ... x[ngp-1] are x1 ... xN) 
//    ekin    : DVR of the kinetic energy operator   
//
//
void Tdiag(int ngp, double mu, double x0, double xNp1, double *x, double *ekin, int ldt, int verbose)
{
  double DeltaX = (xNp1 - x0) / (double)(ngp+1);   // L = x(N+1) - x0

 // verbose = 4;
  if (verbose > 0) {
    printf("FFT  with %d grid points\n", ngp);
    printf("  mass: %f\n", mu);
    printf("  box walls at %10.6f  -  %10.6f,  L = %10.6f\n", x0, xNp1, xNp1-x0);
    printf("  grid range   %10.6f  -  %10.6f\n", x0+DeltaX, xNp1-DeltaX);
  }

  
  //  note that x0 and xN+1 are not grid points
  //
  int j;
  for (j = 0; j < ngp; ++j) {  
    x[j] = x0 + DeltaX * (double)(j+1);
   // cout<<" x["<<j<<"]= "<<x[j]<<endl;
  }
  if (verbose > 1) {
    printf("grid points:\n");
    for (int j = 0; j < ngp; ++j)
      printf("%3i  %12.6f\n", j+1, x[j]);
  }

  
  double f1 = 0.5 * PI*PI / (mu * SQR(DeltaX));
  double f2 = PI / (double)(ngp);
  double f3 = SQR(ngp+1);
  int a, b, ic;

  //
  //  diagonal terms
   /*************new part*************/
   double DeltaK = 2.0*f2/DeltaX;    
  // cout<<"DeltaK = "<<DeltaK<<endl;
    for( a = 0; a < ngp; a++){
      int kpoint  = a;
     //  kpoint=-ngp/2+a;
      if(kpoint >= ngp/2) kpoint= kpoint - ngp; 
      // k[x]        = dk*kpoint;
      ekin[a]      = pow(DeltaK*kpoint,2)/(2.0*mu);
    // cout<<"igp and kpoint ekin="<<a<<"  "<<kpoint<<"  "<<ekin[a]<<endl;
    }
    /************************************/ 
   
}
