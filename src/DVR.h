#ifndef PISCES_DVR_H_
#define PISCES_DVR_H_

#include <iostream>
#include "vecdefs.h"

//******header files for FFT
// header files for FFTW3
//#include "/opt/sam/fftw/3.3.3/intel/openmp/fast/include/fftw3.h"
#include "fftw3.h"
//#include "/ihome/crc/install/intel-2017.1.132/intel-mpi-2017.1.132/fftw/3.3.5/include/fftw3-mpi.h"
//#include "fftw3-mpi.h"

//intel mkl wrapper for fftw3
//#include "fftw3_mkl.h"

#include <complex>
typedef std::complex<double> Complex;
//******************************

/** \brief Provides functions to construct a %DVR wavefunction for the excess electron.
    
    A \em %DVR means the wavefunction is on a grid, however, the grid is based on underlying 
    basis functions, here either Sine=particle-in-a-box or Harmonic Oscillator
*/


class DVR 
{
public:
   /// Initializes data members to zero
   DVR() 
      : verbose(0)
      , no_dim(MAXDIM)
      , tformat(0) 
      , nconverged(0)
      , nwavefn(0)
      , StepSize(MAXDIM)
   {}

   /// Deallocates work arrays
   ~DVR();


   void SetVerbose(int v);


   /** \brief Set up grid.

   \param nx  Number of grid points in x-dir
   \param ny  Number of grid points in y-dir
   \param nz  Number of grid points in z-dir
   \param type 1 for HO, 2 for Sine DVR
   \param gpara Three frequencies (HO) or three lengths (Sine). Both in au
   */
   void SetupDVR(const int *nx, int type, int Sampling, const double *gpara, int verbose);
   void SetupDVR2(const int *nx, int type, int Sampling, const double *gpara, int verbose);
   
   
   /// \brief Evaluates the potential at every grid point
   ///
   /// This is the DVR approximation: the V-operator is diagonal.  
   void ComputePotential(class Potential &V);

//   void ComputeGradient(class Potential &V, int nSites, double *Gradient, double *DmuByDx, double *DmuByDy, double *DmuByDz);
//   void ComputeGradient(class Potential &V, int nSites, double *Gradient, double *dEfield, double *PolGrad, class WaterCluster &WaterN);
   void ComputeGradient(class Potential &V, int nSites, double *Gradient, double *dT_x , double *dT_y , double *dT_z, double *PolGrad, class WaterCluster &WaterN , double *dEfield);



   /** \brief Calls an iterative Eigen-solver (Lanczos-Arnoldi or Davidson) to compute the energy and wavefunction of the excess electron
   
   Both are subspace methods. Call DiagonalizeSetup() first. 

   \param SVFlag    start vector flag
   \li 1 : lowest particle-in-the-box startvector
   \li 0 : the start vector is read from a smaller grid (StartVector.cub)
   \li 2 : to use the old wavefunction 
   */ 
   int Diagonalize(int SVFlag, double *ev);

   /** \brief Set parameters for the Eigen-solver
   
   \param nEV  no of Eigen-pairs (Davidson can do only 1 so far)
   \param DiagFlag 1=Arnoldi 2=Davidson 3=zeroth-order Jacobi-Davidson
   \param nMaxSub  maximal no of vectors in subspace expansion
   \param nMaxIter  no of macro iterations (a micro iteration is one matrix-times-vector,
                    a macro iteration is nMaxSub matrix-times-vector operations)
   \param pTol      convergence tolerance is 10^{-pTol} 
   */
   void DiagonalizeSetup(int nEV, int DiagFlag, int nMaxSub, int nMaxIter, int pTol);

   /// Writes a gaussian-type cubefile that can be read by Avogadro/
   void WriteCubeFile(int iwf, const char *fname, int nAtoms, const int *Z, const double *position, int cubeflag = 0);

   /// put the DVR wavefunctions into cube 
   void GetWaveFnCube(int iwf, double *cube);

   /// Writes 2D grids 2D grids for plotting the wavefunction
   void WriteCuts(void);
   void WriteOneDCuts(void);
   
   void ExpectationValues(int verbose=0);
   void EnergyPartitioning(class Potential &V) ;
private: // METHODS

   void ComputeGridParameters();
   void ComputeGridPointsAndKineticEnergy();
   //set up the FFT parameters
   void FFTSetup();
   /// Can be used as a start vector
   void ParticleInAnDBoxWf(double *wf);  
   void ExtendWfDouble(double *wf);  
   void ExtendWfThree(double *wf);  
   double InterpolVdouble(double *TempF, int px, int py, int pz, int *Pre1db);
   double InterpolVtriple(double *TempF, int px, int py, int pz, int *Pre1db);


   /// \name Diagonalizer functions
   //@{
   void MatrixTimesVector(const double *x, double *y);
   //peforms the FFT part instead of matrix times vector 
    void VectorFFT_old(const double *x, double *y);
   int larnoldi(int ng, int nev, int maxsub, int maxiter, int ptol, double *ev);
   int davdriver(int ng, int nstates, int maxsub, int maxiter, int ptol, int jdflag, double *ev);
   void ComputeDiagonal(double *diag);
   // for debugging a full diagonalization 
   void build_h(double *hmat);
   void fulldiag(double *hop);
   //@}

private: // MEMBERS 
   class counter {
   public:
     counter () : c(0),verbose(0) {}
     ~counter () { if (this->verbose > 0) std::cout << "MatrixTimesVector-count " << c << std::endl; }
      int c;
      int verbose;
   };
   counter* m_pkc;
   static const int MAXDIM = 3;
   int verbose;
   int dvrtype;          ///< 1 : Harmonic oscillator; !=1 Sine
   int sampling;         // evaluate V at gridpoints, or double, or triple density
   double gridpara[MAXDIM];///< gridparameter: for type 1: omega; for type 2: length
   int no_dim;           ///< no of dimensions (3 here, but can be arbitrary but many functions make sense only for 3) 
//   int max1db;           ///< maximum of grid points in any dimension
   int max1db[3];           ///< maximum of grid points in any dimension
//   int Pre1db;           ///< previous maximum of grid points in any dimension
   int Pre1db[3];           ///< previous maximum of grid points in any dimension
   int n_1dbas[MAXDIM];  ///< no of grid points for each dimension
   int incv[MAXDIM];     ///< stride for the grid
   int ngp;              ///< no of grid points
   int tformat;          ///< format of T (triangle=0 or full matrix=1)
   int nconverged;       ///< no of converged eigenvalues and eigenfunctions (from ARPACK or Davidson)
   int nwavefn;          ///< current number of wavefunctions allocated
   int DiagCount;

   ///\name diagonalization information
   //@{
   int maxSub;   ///< max no of micro-iterations
   int maxIter;  ///< max no of macro-iterations 
   int ptol;     ///< tolerance = 10^-pTol 
   int nStates;  
   int diagFlag; ///< 1=Lanczos  2=Davidson  3= 0th-order Jacobi-Davidson
   //@}

   int Idual;  // check the nubmer of grid is odd or even
   double* e_kin[MAXDIM];   ///< kinetic energy matrices of coordinate i
   double* dvr_rep[MAXDIM]; ///< matrix to go from DVR to FBS representation              
   dVec Vec_x_dvr;
   double*  x_dvr;   ///< list with gridpoints in each dimension : x_dvr[max_1db * no_dim]
   dVec Vec_v_diag;  
   double*  v_diag;  ///< potential energy at the grid points v_diag[product of n_1dbas]
   dVec Vec_v_diag_pc;  
   double*  v_diag_pc;  
   dVec Vec_v_diag_ind;  
   double*  v_diag_ind;  
   dVec Vec_v_diag_rep;  
   double*  v_diag_rep;  
   dVec Vec_v_diag_pol;  
   double*  v_diag_pol;  
   //double*  wavefn;
   dVec wavefn;  ///< value of the wavefunction at the grid points (times volume element)
   dVec CoarseWf;  ///< value of the wavefunction at the sparse grid points (times volume element)
   dVec StepSize;

//   iVec select;       // the bloody, allegedly not referenced array
//   dVec v;             // Lanczos basis
//   dVec workd;         // Lanczos vectors
//   dVec workl;    // work space
//   dVec resid;    // residual vector
   int iparam[11], ipntr[14];

   int SaveIparam[11];
// SaveIpntr[14], SaveNev;

   dVec SaveV;             // Lanczos basis
   dVec SaveWorkd;         // Lanczos vectors
   dVec SaveWorkl;    // work space
   dVec SaveResid;    // residual vector

   
   //Variable for the fourier transformation of x and KE operator  
   int threads_ok;
   double* KE_diag;   ///< final diagonal kinetic energy matrices coordinate x,y,z

  Complex *phi_x;     //Copies vector x to new complex array: x can be directly fourier transformed
  Complex *phi_k;    //forward fourier transform of vector x
  Complex  *KE_phi_k; //KE matrix multiplied by FFT(x):T*FFT(x)
  Complex  *KE_phi_x; //  Backward fourier transform of T*FFT(x)
 //  fftw_complex *phi_x;
 //  fftw_complex  *phi_k;
//   fftw_complex *KE_phi_x;
 //  fftw_complex *KE_phi_k;

   fftw_plan plan_forward;
   fftw_plan plan_backward;
   int xa1,ya2,za3;
};


#endif // PISCES_DVR_H_
