#include <iostream>
#include <math.h>
using namespace std;

///
///  this takes a wavefunction on a small grid n[3], xmin[3], dx[3], wf[] 
///  and creates a rough image on a large grid (capitals) to be used 
///  as a startvector for a Lanczos or Davidson
///
///  this can easily be extended to find the nearest 8 indices (ixl, ixu, iyl, ...)
///  and compute a weighted average of those 8 points
///  however, for a start vector a bit of random noise seems more healthy
///
///
void Small2Large(int *n, double *xmin, double *dx, double *wf,
		 int *N, double *XMIN, double *DX, double *WF)
{

  int verbose = 1;
  int CubeFlag = 1;   // 1 wf is from a cube; 0: wf is a DVR wavefunction (0 not implemented yet)
  int InterpolationFlag = 1; // 0: nearest   1: eight weighted points 


  if (verbose > 0) {
    cout << "Small grid:\n"
	 << n[0] << " x pts from " << xmin[0] << " to " << xmin[0]+(n[0]-1)*dx[0] << " steps of " << dx[0] << "\n" 
	 << n[1] << " y pts from " << xmin[1] << " to " << xmin[1]+(n[1]-1)*dx[1] << " steps of " << dx[1] << "\n" 
	 << n[2] << " z pts from " << xmin[2] << " to " << xmin[2]+(n[2]-1)*dx[2] << " steps of " << dx[2] << "\n" 
	 << "Large grid:\n"
	 << N[0] << " x pts from " << XMIN[0] << " to " << XMIN[0]+(N[0]-1)*DX[0] << " steps of " << DX[0] << "\n" 
	 << N[1] << " y pts from " << XMIN[1] << " to " << XMIN[1]+(N[1]-1)*DX[1] << " steps of " << DX[1] << "\n" 
	 << N[2] << " z pts from " << XMIN[2] << " to " << XMIN[2]+(N[2]-1)*DX[2] << " steps of " << DX[2] << "\n"; 
    cout << "1st elements of wf:" << wf[0] << " " << wf[1] << "\n";

  }



  int ix = 0, iy = 0, iz = 0;
  double xmax[3] = {0,0,0};
  for (int k = 0; k < 3; ++k)
    xmax[k] = xmin[k] + (n[k]-1)*dx[k];
  int xstride = 1;
  int ystride = 1;
  int zstride = 1;
  if (CubeFlag) {
    // small grid has been read from a cube:
    // z is fastest, zstride = 1
    ystride = n[2];
    xstride = n[1] * n[2];
  }

    
  // loop over the LARGE grid
  // internally the first dimension is fastest (has stride 1) 
  //
  double *pWF = WF;
  double nrm = 0;
  double value = 0;


  if (InterpolationFlag == 0)
    {
      
      for (int IZ = 0; IZ < N[2]; ++IZ)
	{
	  double Z = XMIN[2] + IZ * DX[2];
	  if (Z < xmin[2])
	    iz = 0;
	  else if (Z > xmax[2])
	    iz = n[2]-1;
	  else
	    iz = int(0.49 + (Z-xmin[2])/dx[2]);
	  
	  for (int IY = 0; IY < N[1]; ++IY) 
	    {
	      
	      double Y = XMIN[1] + IY * DX[1];
	      if (Y < xmin[1])
		iy = 0;
	      else if (Y > xmax[1])
		iy = n[1]-1;
	      else
		iy = int(0.49 + (Y-xmin[1])/dx[1]);
	      
	      for (int IX = 0; IX < N[0]; ++IX)  
		{
		  double X = XMIN[0] + IX * DX[0];
		  if (X < xmin[0])
		    ix = 0;
		  else if (X > xmax[0])
		    ix = n[0]-1;
		  else
		    ix = int(0.49 + (X-xmin[0])/dx[0]);
		  
		  value = wf[ix*xstride + iy*ystride + iz*zstride];
		  
		  
		  //debug
		  // double x = xmin[0] + ix * dx[0];
		  // double y = xmin[1] + iy * dx[1];
		  // double z = xmin[2] + iz * dx[2];
		  // cout << "Approximating Psi at " 
		  //      << "( " << X << ", " << Y << ", " << Z << ") with psi at " 
		  //      << "( " << x << ", " << y << ", " << z << ") = " << value << "\n"; 
		  
		  *pWF = value;
		  nrm += value * value;
		  pWF += 1;
		}
	    }
	}
    }
  else  // trilinear interpolation from eight corners of the surrounding cube
    {
      int ixu, iyu, izu;
      double wx, wy, wz;

      for (int IZ = 0; IZ < N[2]; ++IZ)
	{
	  double Z = XMIN[2] + IZ * DX[2];
	  iz = int((Z-xmin[2])/dx[2]);
	  if (iz < 0) 
	    {iz = 0; izu = 1;}
	  else if (iz > n[2]-2)
	    {iz = n[2]-2; izu = n[2]-1;}
	  else
	    izu = iz + 1;
	  wz = (Z - (xmin[2] + iz*dx[2])) / dx[2];


	  for (int IY = 0; IY < N[1]; ++IY) 
	    {
	      double Y = XMIN[1] + IY * DX[1];
	      iy = int((Y-xmin[1])/dx[1]);
	      if (iy < 0) 
		{iy = 0; iyu = 1;}
	      else if (iy > n[1]-2)
		{iy = n[1]-2; iyu = n[1]-1;}
	      else
		iyu = iy + 1;
	      wy = (Y - (xmin[1] + iy*dx[1])) / dx[1];
	   
	      for (int IX = 0; IX < N[0]; ++IX)  
		{
		  double X = XMIN[0] + IX * DX[0];
		  ix = int((X-xmin[0])/dx[0]);
		  if (ix < 0) 
		    {ix = 0; ixu = 1;}
		  else if (ix > n[0]-2)
		    {ix = n[0]-2; ixu = n[0]-1;}
		  else
		    ixu = ix + 1;
		  wx = (X - (xmin[0] + ix*dx[0])) / dx[0];


		  // eight points:
		  double sumdxyz = 0;
		  double weight = 0;
		  weight = (1-wx) * (1-wy) * (1-wz);
		  sumdxyz += weight;
		  value = wf[ix*xstride + iy*ystride + iz*zstride] * weight;
		  
		  weight = wx * (1-wy) * (1-wz);
		  sumdxyz += weight;
		  value += wf[ixu*xstride + iy*ystride + iz*zstride] * weight;

		  weight = (1-wx) * wy * (1-wz);
		  sumdxyz += weight;
		  value += wf[ix*xstride + iyu*ystride + iz*zstride] * weight;

		  weight = (1-wx) * (1-wy) * wz;
		  sumdxyz += weight;
		  value += wf[ix*xstride + iy*ystride + izu*zstride] * weight;

		  weight = wx * wy * (1-wz);
		  sumdxyz += weight;
		  value += wf[ixu*xstride + iyu*ystride + iz*zstride] * weight;

		  weight = wx * (1-wy) * wz;
		  sumdxyz += weight;
		  value += wf[ixu*xstride + iy*ystride + izu*zstride] * weight;

		  weight = (1-wx) * wy * wz;
		  sumdxyz += weight;
		  value += wf[ix*xstride + iyu*ystride + izu*zstride] * weight;

		  weight = wx * wy * wz;
		  sumdxyz += weight;
		  value += wf[ixu*xstride + iyu*ystride + izu*zstride] * weight;
		  
		  value /= sumdxyz;

		  *pWF = value;
		  nrm += value * value;
		  pWF += 1;
		}
	    }
	}
    }



  // normalize as DVR wavefunction

  pWF = WF;
  nrm = 1.0 / sqrt(nrm);

  cout << "DVR normalization factor is " << nrm << "\n";

  for (int IZ = 0; IZ < N[2]; ++IZ) {
    for (int IY = 0; IY < N[1]; ++IY) {  
      for (int IX = 0; IX < N[0]; ++IX) {  
	*pWF *= nrm;
	pWF += 1;
      }
    }
  }

}
