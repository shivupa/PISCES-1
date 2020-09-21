//  
//
//  this routine performs  y = H * x
//  where H is a DVR of an n-D Hamiltonian and x and y are wavefunctions on the grid
//
//
#include "DVR.h"
#include "lapackblas.h"



// small inline utility to convert from 3D subscript into Fortran style (column-major) indices
inline int sub2ind( int i, int j, int k, const int n[] )
{
   enum {X,Y,Z};
   return i + n[X]* (j + n[Y]*(k));
}


void DVR::MatrixTimesVector(const double *x, double *y)
{
   m_pkc->c++;

   //
   //  y = V*x  where V is a diagonal matrix (this initializes y at the same time)
   //
#pragma omp parallel 
{
   #pragma omp for
   for (int igr = 0; igr < ngp; igr++)
      y[igr] = v_diag[igr] * x[igr];

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

   // X dimension
   #pragma omp barrier
   #pragma omp for
   for (int k=0; k<n[Z]; ++k) {
      for (int j=0; j<n[Y]; ++j)
      {
         int offset = sub2ind(0,j,k,n);
         dspmv("U", n[X], 1, e_kin[X], &x[offset], stride[X], 1, &y[offset], stride[X]); 
      }
   }
   #pragma omp barrier
   
   // Y dimension
   #pragma omp for
   for (int k=0; k<n[Z]; ++k) {
      for (int i=0; i<n[X]; ++i)
      {
         int offset = sub2ind(i,0,k,n);
         dspmv("U", n[Y], 1, e_kin[Y], &x[offset], stride[Y], 1, &y[offset], stride[Y]); 
      }
   }

   
   // Z dimension
   #pragma omp barrier
   #pragma omp for
   for (int j=0; j<n[Y]; ++j) {
      for (int i=0; i<n[X]; ++i)
      {
         int offset = sub2ind(i,j,0,n);
         dspmv("U", n[Z], 1, e_kin[Z], &x[offset], stride[Z], 1, &y[offset], stride[Z]); 
      }
   }

} // pragma omp parallel
}
