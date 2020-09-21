//  reading Gaussian-cube-file format

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cctype>
#include <iostream>
#include <fstream>

using namespace std;


//  const double Angs2Bohr = 1.889725989;


void ReadCubeFile(int *npts, double *x0, double *dx, double **cube)
{

  cout << "\nReading a start vector from cube file\n";

  int i;
  int verbose = 1;
  const char *fngrid = "StartVector.cub";
  ifstream input;
  input.open(fngrid);


  const int LineLength = 132;
  char buffer[LineLength+1];

  // two comment lines, after that formatted input
  input.getline(buffer, LineLength);
  input.getline(buffer, LineLength);

  
  // read no of atoms and grid parameters
  int nnuc;
  int nx, ny, nz;
  double xx, yy, zz;

  input.getline(buffer, LineLength);
  sscanf(buffer, "%i %lf %lf %lf", &nnuc, x0, x0+1, x0+2);
  //input >> nnuc >> x0[0] >> x0[1] >> x0[2];

  input.getline(buffer, LineLength);
  sscanf(buffer, "%i %lf %lf %lf", &nx, &(dx[0]), &yy, &zz);
  //  input >> nx >> dx[0] >> yy >> zz;
  if (yy != 0 || zz != 0) {cout << "Only cartesian grids please:"<< buffer <<"\n"; exit(1);}

  input.getline(buffer, LineLength);
  sscanf(buffer, "%i %lf %lf %lf", &ny, &xx, &(dx[1]), &zz);
  //  input >> ny >> xx >> dx[1] >> zz;
  if (xx != 0 || zz != 0) {cout << "Only cartesian grids please:"<< buffer <<"\n"; exit(1);}

  input.getline(buffer, LineLength);
  sscanf(buffer, "%i %lf %lf %lf", &nz, &xx, &yy, &(dx[2]));
  //  input >> nz >> xx >> yy >> dx[2];
  if (xx != 0 || yy != 0) {cout << "Only cartesian grids please:"<< buffer <<"\n"; exit(1);}
  npts[0] = nx;
  npts[1] = ny;
  npts[2] = nz;

  int ntotal = nx * ny * nz;
  *cube = new double[ntotal];
  
  if (verbose > 0) {
    cout << nx << " x pts from " << x0[0] << " to " << x0[0]+(nx-1)*dx[0] << " steplength = " << dx[0] << endl;
    cout << ny << " y pts from " << x0[1] << " to " << x0[1]+(ny-1)*dx[1] << " steplength = " << dx[1] << endl;
    cout << nz << " z pts from " << x0[2] << " to " << x0[2]+(nz-1)*dx[2] << " steplength = " << dx[2] << endl;
  }

  // read over the lines with the nuclei
  for (i = 0; i < nnuc; ++i) {
    input.getline(buffer, LineLength);
  }

  // 
  //  read the grid data and do your sums  
  //
  double density = 0; 

  int nvaluesread = 0;
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
	// get next density
	input >> density;
	(*cube)[nvaluesread] = density;
	if (nvaluesread < 3) {
	  cout << nvaluesread << "-th value read from cube file: " << density << "\n";  
	}
	nvaluesread += 1;
      }
    }
  }
  input.close();
  printf("Last value read is %e (no %i)\n", density, nvaluesread);
  if (nvaluesread != ntotal)
    cout << "But there should be " << ntotal << endl;
  
}


