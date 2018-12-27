#include <iomanip>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cassert>

#ifdef PISCES_USE_MPI
#include <mpi.h>
#endif

#include <mpi.h>
#include <fftw3-mpi.h>

#include "timer.hpp"
#include "tsin.h"
#include "constants.h"
#include "vecdefs.h"
#include "Parameters.h"
#include "GetInput.h"

// DPP incudes
#include "GTO.h"
#include "MO.h"
#include "AtomCenter.h"
#include "Water.h"
#include "DPP.h"
#include "Molecule.h"


// excess electron includes
#include "Potential.h"
#include "DVR.h"
#include "ClusterAnion.h"
// Optimization
//#include "optimize.h"

//plotting
#include "polplot.h"

// NaCl includes
#include "NaCl.h"
#include "C60.h"
#include "Model_pot.h"

// potfit includes
#include "potfit.h"

using namespace std;

///////////////////////////////////////////////////////////////////////
//
// Entry point to pisces
//
//
int main(int argc, char* argv[])
{
//  MPI_Init( &argc, &argv );
  int provided, threads_ok;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  threads_ok=provided>=MPI_THREAD_FUNNELED;
  cout<<"threads_ok = "<<threads_ok<<endl;

  int rank, size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Status status;

  if(rank == 0)  cout << "\nPittsburgh Infra Structure for Clusters with excess ElectronS\n";
   progress_timer time("Total");

  if(rank ==0) cout << "MPI: num of MPI ranks = " << size << endl;
   // mpiworkshare();



#ifdef _OPENMP
 if(rank == 0)  cout << "OpenMP: num_threads = " << omp_get_max_threads() << endl;
#endif

  if (argc != 2)  {
    if(rank == 0) cout << "usage: pisces input-file\n";
    return EXIT_FAILURE;
  }

  // read input file and store it line-by-line in TSIN Input
  TSIN Input;
  if(rank == 0) cout << "\nThis is the Input read from " << argv[1] << ":\n";
  Input.ReadFromFile(argv[1], 5);

  // Parse all input groups, and put parameters into InP
  Parameters InP;
  GetInputParameters(InP, Input);
  


  // If dealing with water clusters, get the coordinates of the water monomers 
  dVec WaterCoor;
  if (InP.runtype < 100) { // above is NaCl clusters
    int nW = ReadWatersGroup(Input, WaterCoor);
    if (InP.nWater == 0)
      InP.nWater = nW;
    else if (nW != InP.nWater && InP.nWater > 0) 
       throw std::invalid_argument("inconsistent no of waters in Waters and WaterModel groups");      
  }

  //if dealing with the electron
//  ClusterAnion Wn_test;
  // Can we run it?
  switch (InP.runtype) 
    {
    case 1: 
    if(InP.runtype == 1) {  ClusterAnion Wnn; 
    Wnn.SetUpClusterAnion(WaterCoor, InP);
    Wnn.SinglePoint(WaterCoor);}
      break;
    case 2:
     Optimizer(WaterCoor, InP);
      break;
    case 42:
      potfit(WaterCoor, InP);
      break;
//vkv
    case 60:
      polplot(WaterCoor, InP);
      break;
//vkv
    case 101:
      NaCl_sp(Input);
      break;
    case 102:
      C60_sp(Input);
      break;
    case 103:
      SphericalPotential(Input);
      break;
    default:
       std::ostringstream error;
       error << "unknown runtype: " << InP.runtype;
       throw std::invalid_argument( error.str() );
    }
  
  return EXIT_SUCCESS;
}

//catch( const std::exception& e ) {
//   std::cerr << "error: " << e.what() << std::endl;
//   return EXIT_FAILURE;
//}

