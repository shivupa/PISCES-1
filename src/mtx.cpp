//  
//
//  this routine performs  y = H * x
//  where H is a DVR of an n-D Hamiltonian and x and y are wavefunctions on the grid
//
//
#include "DVR.h"
#include "lapackblas.h"

#include "timer.hpp"



// small inline utility to convert from 3D subscript into Fortran style (column-major) indices
inline int sub2ind( int i, int j, int k, const int n[] )
{
   enum {X,Y,Z};
   return i + n[X]* (j + n[Y]*(k));
}


void DVR::MatrixTimesVector(const double *x, double *y)
{

 progress_timer t("MatrixTimesVector", verbose);
   m_pkc->c++;

   //
   //  y = V*x  where V is a diagonal matrix (this initializes y at the same time)
   //
#pragma omp parallel 
{
 //  cout<<"ngp="<<ngp<<endl;
   #pragma omp for
   for (int igr = 0; igr < ngp; igr++) {
      y[igr] = v_diag[igr] * x[igr];
    // cout<<"x["<<igr<<"]="<<x[igr]<<endl;
   }
   //
   // x and y is [n0][n1][n2]
   // e0 is [n0][n0]
   // e1 is [n1][n1]
   // e2 is [n2][n2]
   // e0,e1,e2 are symmetric
   //
   enum{X,Y,Z};
   // alias (`n' looks cleaner than `n_1dbas`)
   const int (&n) [MAXDIM] = n_1dbas;
   const int stride[] = {1,n[X],n[X]*n[Y]};

   double done = 1.0;
   double dzero = 0.0;
   int one = 1;



   // cout<<"n[X]="<<n[X]<<" n[Y]="<<n[Y]<<" n[z]="<<n[Z]<<endl;
   // X dimension
   #pragma omp barrier
   #pragma omp for
   for (int k=0; k<n[Z]; ++k) {
      for (int j=0; j<n[Y]; ++j)
      {
         int offset = sub2ind(0,j,k,n);
         //cout<<"1X: y["<<offset<<"]= "<<y[offset]<<endl;
         //cout<<"0jk"<<" 0 "<<j<<" "<<k<<endl;
         dspmv("U", n[X], 1, e_kin[X], &x[offset], stride[X], 1, &y[offset], stride[X]); 
         //dsymv("U", &n[X], &done, e_kin[X], &n[X], &x[offset],&one, &done, &y[offset], &one); 
         // cout<<"stride[X]="<<stride[X]<<" *e_kin[X]="<<*e_kin[X]<<endl;
         //cout<<"X="<<X<<" offset="<<offset<<endl;
         //cout<<"2X: y["<<offset<<"]= "<<y[offset]<<endl;
         //cout<<"e_kin["<<X<<"]= "<<e_kin[X]<<endl;
      }
   }
   #pragma omp barrier
   
   // Y dimension
   #pragma omp for
   for (int k=0; k<n[Z]; ++k) {
      for (int i=0; i<n[X]; ++i)
      {
         int offset = sub2ind(i,0,k,n);
        // cout<<"i0k"<<i<<" 0 "<<k<<endl;
        //cout<<"1Y: y["<<offset<<"]= "<<y[offset]<<endl;
         dspmv("U", n[Y], 1, e_kin[Y], &x[offset], stride[Y], 1, &y[offset], stride[Y]); 
        //dsymv("U", &n[Y], &done, e_kin[Y],&n[Y], &x[offset],&one, &done, &y[offset], &one); 
        //cout<<"2Y: y["<<offset<<"]= "<<y[offset]<<endl;
      }
   }

   
   // Z dimension
   #pragma omp barrier
   #pragma omp for
   for (int j=0; j<n[Y]; ++j) {
      for (int i=0; i<n[X]; ++i)
      {
         int offset = sub2ind(i,j,0,n);
         //cout<<"ij0"<<i<<" "<<j<<" 0 "<<endl;
         dspmv("U", n[Z], 1, e_kin[Z], &x[offset], stride[Z], 1, &y[offset], stride[Z]); 
         //dsymv("U", &n[Z], &done, e_kin[Z],&n[Z], &x[offset],&one, &done, &y[offset], &one); 
        //cout<<"Z: y["<<offset<<"]= "<<y[offset]<<endl;
      }
   }



} // pragma omp parallel
}
