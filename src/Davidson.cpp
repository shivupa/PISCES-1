///////////////////////////////////////////////////////////////////////////////
///
///   Davidson 
///  
///   multi-root, reverse interface
///
///   for the notation of the variables:
///
///   we try to compute eigenpairs of the matrix \f$H\f$
///
///   a basis set \f$B\f$ (column vectors of length ndim) is iteratively build up 
///   \f$Z = HB\f$ are the vectors produced by applying \f$H\f$ to the basis vectors
///   \f$S = BtHB = BtZ \f$ is the representation of \f$H\f$ in the subspace spanned by \f$B\f$
///   \f$V\f$ are the eigenvectors of \f$S\f$
///
///   then the eigenvalues of \f$S\f$ are approximations for the eigenvectors of \f$H\f$
///   and \f$BV\f$ are the associated approximations for the eigenvectors of \f$H\f$
///
///////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "lapackblas.h"
#include "vecdefs.h"

//
//  functions used by the Davidson and defined later in this file
//
void ComputeS(int ndim, int nsubsp, int nadd, int maxsubsp, double *B, double *Z, double *S, int verbose);
void OrthoVecOnB(int ndim, int nbas, double *vec, double *B, int verbose);
void GramSchmidt(int ndim, int nvec, double *vec, int verbose);
void DavidsonJacobiCorrectionVector(int ndim, double lambda, double *res_vec, double *jd_vec,
				    double *ritz_vec, double *diagH);
void DavidsonCorrectionVector(int ndim, double lambda, double *res_vec, double *diagH);

using namespace std;


///////////////////////////////////////////////////////////////////////////////
///
///   compute size of work array needed 
///   see Davidson code comments where davwork assigments are made 
///
///////////////////////////////////////////////////////////////////////////////
int DavidsonWorkSize(int ndim, int maxsub, int nroots, int corrflag)
{
  if (corrflag == 2)
    return 2*maxsub*maxsub + 11*maxsub + (2*nroots+1)*ndim;  // one additional vector for Jacobi correction
  else
    return 2*maxsub*maxsub + 11*maxsub + 2*nroots*ndim;
}


///////////////////////////////////////////////////////////////////////////////
///
///  wasteful Davidson: still experimental stage, I want to keep it all,
///  because what used to work great for CI matrices sucks for grids 
///
///////////////////////////////////////////////////////////////////////////////
int Davidson(int ndim,        /// dimension of H, length of vectors in B and Z
	     int maxsub,      /// max subspace dimension
	     int nroots,      /// no of roots of H to be found
	     int maxmacro,    /// max no of macro iterations
	     int tol,         /// norm of ritz vectors must be below 10^-tol
	     int corrflag,    /// correction vector flag (see below)
	     int verbose,     /// verbocity level
             double *evals,   /// returns converged eigenvalues
	     int &nConv,      /// returns number of converged eigenpairs
	     double *B,       /// space for basis set build by Davidson iteration
	     double *Z,       /// space for H*B vectors
	     double *diag,    /// diag(H) for Davidson correction vectors
	     double *davwork, /// work space for subspace arrays and ritz vectors 
	     int *inout)      /// communication codes with the reverse interface
///////////////////////////////////////////////////////////////////////////////
///
///  startvectors are on input in B
///  on output converged ritzvectors are put into B
///
///
///   corrflag = 0  the residual vector is added to B: effectively Lanczos-Arnoldi
///            = 1  standard Davidson correction vector
///            = 2  zeroth-order Jacobi-Davidson correction
///
///      more expensive correction vectors usually mean fewer matrix-times-vector calls
///      where the best tradeoff is case dependent, and a matter of experimentation
///      do not expect the same miracles the Davidson does for GTO sets for grids 
///
///   davwork will hold the following arrays
///
///    S = new double [maxsubspace * maxsubspace];
///    V = new double [maxsubspace * maxsubspace];
///    Work = new double[lwork];   lwork = 10*maxsub
///    sse = new double[maxsubspace];
///    ritzvec  = new double[ndim];
///    if (jdflag == 1) jdvec = new double[ndim];
///
///    size of davwork is 2M^2 + 11M + N (+N for jdflag) 
///    this is computed in DavidsonWorkSize()
///
///    internal status
///      0 : first call, initialize dimensions, set pointers, normalize startvectors, and switch status to 1
///      1 : standard
///
///   return codes:   0: done; converged (eval and B return the converged eigenpairs) or maxmacro exceeded
///                            the number of converged vectors is returned in inout[0]
///                            if maxmacro is exceeded a list of current convergence infromation is printed
///                   1: call Matrix-times-vector, i.e., Z[..] := H * B[..]  
///
///
///////////////////////////////////////////////////////////////////////////////
{

  static int mystatus = 0;
  static double thresh = std::pow(10.,-tol);

  static double *rcvec, *ssvec;  /// these are truly pointers to existing blocks
  static int lwork = 0;

  static int jmacro = 0;     /// no of macro iterations done
  static int nsubspace = 0;  /// current number of subspace vectors
  static int nnewBs = 0;     /// number of vectors selected for the next iterations
  static int nnewZs = 0;     /// number of vectors just procecessed by mtx
  static int newBmax = 1;    /// use nroots for GTO tasks (startspace is meaningful ); 
                             /// use 1 for grids (startspace could as well be random)
  
  static int one = 1;        /// arguments for 
  static double done = 1.0;  /// BLAS and LAPACK 
  static double dzro = 0.0;  /// functions

  static dVec residuals;     /// residuals |r| = norm of residual vectors of each root
  static dVec conv_weights;  /// weight of each subspace vector with already converged space
  static iVec ccode;         /// ritz-vector is converged (0),  newly converged (1), or unconverged (2)
  static iVec Bindex;        /// where in B the ritz vector is put when a marco iteration is started

  // pointers to main internal arrays; these are in davwork
  static double *S = 0;
  static double *V = 0;
  static double *Work = 0;
  static double *sse = 0;
  static double *ritzvecs = 0;
  static double *resvecs = 0;
  static double *jdvec = 0;

  if (mystatus == 0) {
    //
    //  new diagonalization, nothing has been done so far
    //    
    jmacro = 1;  
    nsubspace = 0;
    nnewBs = nroots;
    thresh = pow(10.0, -tol);
    lwork = 10 * maxsub;
    residuals.resize(nroots);
    conv_weights.resize(nroots);
    ccode.resize(nroots);
    Bindex.resize(nroots);
    if (verbose > 0) {
      cout << "Davidson called\n";
      cout << "  Matrix dimension " << ndim << "\n";
      cout << "  Number of roots to be found is " << nroots << "\n";
      cout << "  Maximal subspace size " << maxsub << "\n";
      cout << "  Maximal no of macro iterations " << maxmacro << "\n";
      switch (corrflag)
	{
	case 0:
	  cout << "  Residual vectors are used as correction vectors (effective Lanczos-Arnoldi)\n";
	  break;
	case 1:
	  cout << "  Standard Davidson correction vectors\n";
	  break;
	case 2:
	  cout << "  Simplified Jacobi-Davidson correction vectors\n";
	  break;
	default:
	  cout << "Error in Davidson:  corrflag = " << corrflag << " this should not happen; corrflag set to 1.\n";
	  corrflag = 1;
      }
    }
    //  set pointers to davwork
    S = davwork;                      //new double[maxsub * maxsub];
    V = S + maxsub * maxsub;          //new double[maxsub * maxsub];
    Work = V + maxsub * maxsub;       //new double[lwork];   
    sse = Work + lwork;               //new double[maxsub];
    ritzvecs = sse + maxsub;          //new double[ndim*nroots];
    resvecs = ritzvecs + ndim*nroots; //new double[ndim*nroots]; 
    jdvec = resvecs + ndim*nroots;    //new double[ndim]; only used if jdvec == 1

    //  orthonormalize start vectors in B
    GramSchmidt(ndim, nroots, B, verbose);

    mystatus = 1;
  }
  else if (mystatus != 1) {
    cout << "mystatus = " << mystatus << "This should not happen.\n";
    exit(1);
  }


  //
  //   start of the big while loop
  //

  while (mystatus) {
    //
    //   this is the entry point for the next iteration and for the return from the mtx function
    //   so either return to multiply the new basis vector with H, or pass on into the loop  
    //
    if (nnewBs > 0) {
      if (verbose > 5) {cout << "Davidson calls matrix-times-vector; vectors to multiply = " << nnewBs << "\n"; cout.flush();}
      inout[0] = nsubspace;
      inout[1] = nsubspace;
      inout[2] = nnewBs;
      nnewZs = nnewBs;
      nnewBs = 0;
      return 1; // code to do mtx with inout[2] vectors sittings in B[inout[0]]
    }
    
    //  compute and diagonalize S
    ComputeS(ndim, nsubspace, nnewZs, maxsub, B, Z, S, verbose);
    nsubspace += nnewZs;
    int some_int = maxsub * maxsub; // used for maxsub^2 and for info
    dcopy(&some_int, S, &one, V, &one);
    dsyev("V", "U", &nsubspace, V, &maxsub, sse, Work, &lwork, &some_int);
    if (some_int != 0) 
      { cout << "Davidson diagonalizing S: Lapack error code in DSYEV = " << some_int << endl; exit(1); }

    // compute conv_weights 
    for (int iroot = 0; iroot < nroots; ++iroot) {
      if (nConv > 0) {
	conv_weights[iroot] = dnrm2(&nConv, &V[iroot*maxsub], &one);
	if (verbose > 5)
	  cout << "    E=" << sse[iroot] << " has overlap of " << conv_weights[iroot] << " with converged space.\n"; 
      }
      else
	conv_weights[iroot] = 0;
    }
    
    // build ritz vectors u=Bv and keep for use in Davidson-Jacobi or for starting new marco iteration 
    dgemm("N", "N", ndim, nroots, nsubspace, done, B, ndim, V, maxsub, dzro, ritzvecs, ndim);
    // renormalize ritz-vectors; norms of Bs and Vs are OK, but those of ritz vectors have some noise at the 1e-5 level 
    for (int iroot = 0; iroot < nroots; ++ iroot) {
      double nrm = dnrm2(&ndim, &ritzvecs[iroot*ndim], &one);
      nrm = 1.0 / nrm;
      dscal(&ndim, &nrm, &ritzvecs[iroot*ndim], &one);
    }


    // compute residual vectors r = (H-lambda)Bv = (Z-lambdaB)v and residuals = |r|
    // (replace dgemv and daxpy calls by two dgemm? calls)
    for (int iroot = 0; iroot < nroots; ++iroot) {
      double lambda = sse[iroot];
      ssvec = V + iroot*maxsub;
      rcvec = resvecs + iroot*ndim;
      // step 1: rcvec = Zv = Z * ssvec
      dgemv("N", &ndim, &nsubspace, &done, Z, &ndim, ssvec, &one, &dzro, rcvec, &one);
      // step 2: rcvec = rcvec - lambda ritzvec
      double mlambda = -lambda;
      daxpy(&ndim, &mlambda, &ritzvecs[iroot*ndim], &one, rcvec, &one);
      double curr_res = dnrm2(&ndim, rcvec, &one);
      residuals[iroot] = curr_res;
      if (verbose > 1) {
	if (iroot == 0)
	  printf("%3i   %2i   %16.12f  %10.3e\n", nsubspace, iroot, lambda, curr_res); 
	else
	  printf("      %2i   %16.12f  %10.3e\n", iroot, lambda, curr_res); 
	cout.flush();
      }
    }
    
    // check ritz vectors for convergence (using norm of the residual vector)
    int nnewlyconv = 0;
    double w_converged = 0.95;    
    for (int iroot = 0; iroot < nroots; ++iroot) { 
      if (conv_weights[iroot] > w_converged) {
	ccode[iroot] = 0; // already in converged space, do not correct
      }
      else {
	if (residuals[iroot] < thresh) {
	  ccode[iroot] = 1; // newly converged 
	  Bindex[iroot] = nConv;  // this is its index in B
	  dcopy(&ndim, &ritzvecs[iroot*ndim], &one, &B[nConv*ndim], &one);
	  nConv += 1;
	  nnewlyconv += 1;
	}
	else
	  ccode[iroot] = 2; // unconverged; gets added to B later
      }
    }
    if (nConv == nroots) {
      // Hurrah
      dcopy(&nroots, sse, &one, evals, &one); 
      int all = ndim * nroots;
      dcopy(&all, ritzvecs, &one, B, &one);
      mystatus = 0;
      return 0;
    }

    //
    // new macro-iteration if maxsub exceeded or if at least one vector converged
    // the new macro iteration is started with the rest of the current ritz vectors
    //
    if (nsubspace == maxsub || nnewlyconv > 0) {
      nnewBs = 0;
      // put all unconverged ritz-vectors (nroots-nConv) into B behind the converged ones
      for (int iroot = 0; iroot < nroots; ++iroot) {
	if (ccode[iroot] == 2) {
	  dcopy(&ndim, &ritzvecs[iroot*ndim], &one, &B[(nnewBs+nConv)*ndim], &one);	
	  Bindex[iroot] = nConv+nnewBs;
	  nnewBs += 1;
	}
      }
      if (nnewBs + nConv != nroots) {
	cout << "Warning in Davidson, something is off:\n";
	cout << "nroots=" << nroots << ",  nConv=" << nConv << ", and nnewBs=" << nnewBs << "\n"; 
      }
      // build Z = Z(current)v for the newly converged and unconverged ritz vectors
      // this is build in ritzvecs and then copied as a block into Z
      for (int iroot = 0; iroot < nroots; ++iroot) { 
	if (ccode[iroot] > 0) {
	  int jz = Bindex[iroot];
	  ssvec = V + iroot*maxsub;
	  rcvec = ritzvecs + jz*ndim;  
	  dgemv("N", &ndim, &nsubspace, &done, Z, &ndim, ssvec, &one, &dzro, rcvec, &one);	    
	}
      }
      int all = (nnewBs+nnewlyconv)*ndim;
      dcopy(&all, &ritzvecs[(nConv-nnewlyconv)*ndim], &one, &Z[(nConv-nnewlyconv)*ndim], &one);

      nsubspace = nroots;
      jmacro += 1;
      if (jmacro >= maxmacro) {
	// too many macro iterations have been done
	cout << "\nWarning: The requested number of roots has not been converged after " 
	     << jmacro << " macro iterations.\n";
	for (int iroot = 0; iroot < nConv; ++iroot)
	  printf("%3i  %16.8e  Converged\n", iroot, evals[iroot]);
	for (int iroot = 0; iroot < nroots; ++iroot)
	  printf("%3i  %16.8e  |r|=%10.3e\n", iroot, sse[iroot], residuals[iroot]);
	return 0;
      }
      else {
	if (verbose > 0)
	  cout << "----------- macro iteration " << jmacro << " --------------------\n";
      }
    }

    //
    // add correction vectors to B 
    // (either there was space or a new macro-iteration was started)
    //

    nnewBs = 0;
    for (int iroot = 0; iroot < nroots; ++iroot) {
      if (nsubspace + nnewBs == maxsub)
	break;
      if (ccode[iroot] < 2)
	continue;
      rcvec = resvecs + iroot*ndim;
      if (corrflag == 1)
	DavidsonCorrectionVector(ndim, sse[iroot], rcvec, diag); 
      else if (corrflag == 2)
	DavidsonJacobiCorrectionVector(ndim, sse[iroot], rcvec, jdvec, &ritzvecs[ndim*iroot], diag);
      else 
	; // do not do anything to the residual vector
      // orthogonalise rcvec on the on the vectors already in B

      OrthoVecOnB(ndim, nsubspace+nnewBs, rcvec, B, verbose);
      dcopy(&ndim, rcvec, &one, &B[(nsubspace+nnewBs)*ndim], &one); 
      nnewBs ++;
    }



  } // end of the big while loop
}

///////////////////////////////////////////////////////////////////////////////
///
///   Below is the old code; one needs to compute only the new elements of S, but 
///   as long as B and Z are in memory anyway, recomputing all of S may even be faster
///   and for trying to build a block-algorithm this is a safe fall-back
///
///////////////////////////////////////////////////////////////////////////////
void ComputeS(int ndim, int nsubsp, int nadd, int maxsubsp, double *B, double *Z, double *S, int verbose)
{
  double dzero = 0.0;
  double done = 1.0;
  // S = Zt * B   (see lapackblas.h)

  nsubsp += nadd;

  dgemm("T", "N", nsubsp, nsubsp, ndim, done, Z, ndim, B, ndim, dzero, S, maxsubsp);
  if (verbose > 5) {
    cout << "The S matrix is now:\n";
    for (int i = 0; i < nsubsp; ++i) {
      cout << "S(" << i << "," << i << ") =" << S[i*maxsubsp+i] << "\n";
      if (verbose > 15) {
	for (int j = 0; j < i; ++j) {
	  cout << "S(" << i << "," << j << ") =" << S[i*maxsubsp+j] << " ij-ji=" << S[j*maxsubsp+i]-S[i*maxsubsp+j] <<"\n";
	}
      } 
    }
  }

  //for (int i = 0; i < nsubsp; ++i) 
  //  for (int j = 0; j < i; ++j) 
  //    S[j*maxsubsp+i] = S[i*maxsubsp+j] = 0.5*(S[j*maxsubsp+i]+S[i*maxsubsp+j]);

}



/*    //
    //   compute the new column of S
    //
    if (verbose > 9) {cout << "computing new column of S\n"; cout.flush();}
    if (verbose > 15) cout << "New S elements:";

    nsubpo = nsubspace + 1;
    dgemv("T", &ndim, &nsubpo, &done, B, &ndim, newz, &one, &dzro, S+nsubspace*maxsub, &one);
    for (k = 0; k <= nsubspace; ++k) {
      S[k*maxsub+nsubspace] = S[nsubspace*maxsub+k];
      if (verbose > 15) 
	cout << "  " << S[nsubspace*maxsub+k];
    }
    //  this does it vector by vector:
    //
    //for (k = 0; k <= nsubspace; ++k) {
    //  vec = B + k * ndim;
    //  double s = FORTRANNAME(ddot)(&ndim, vec, &one, newz, &one); 
    //  S[nsubspace*maxsub+k] = S[k*maxsub+nsubspace] = s;
    //  if (verbose > 15) cout << "  " << s;
    //}
    if (verbose > 15) cout << endl;
    nsubspace++;
*/



///////////////////////////////////////////////////////////////////////////////
///
///  orthonormalize the vector vec wrt the basis \f$B\f$
///  it is done twice as the Davidson by construction produces near linear dependent vectors
///
///////////////////////////////////////////////////////////////////////////////
void OrthoVecOnB(int ndim, int nbas, double *vec, double *B, int verbose)
{
  if (nbas == 0) 
    return;
  const int RepeatIt = 2;
  int one = 1;
  if (verbose > 9)
    cout << "  Orthonormalize a new vector " << RepeatIt << " times on " << nbas << " old vectors\n";
  double nrm = 0;
  // first normalize the new vector (vectors in B are assumed to be normalized)
  nrm = 1.0 / dnrm2(&ndim, vec, &one);
  dscal(&ndim, &nrm, vec, &one); 
  // now orthogonalize twice on every vector in B
  for (int i = 0; i < nbas; ++i) {
    double *bv = B + i*ndim;
    for (int doit = 0; doit < RepeatIt; ++doit) {
      double ovl = -ddot(&ndim, vec, &one, bv, &one);
      if (verbose > 10)
        cout << "    overlap " << i << " = " << fabs(ovl) << "\n";
      daxpy(&ndim, &ovl, bv, &one, vec, &one);
      nrm = dnrm2(&ndim, vec, &one);
      if (verbose > 9)
	cout << "    norm-loss-" << doit+1 << "  = " << 1.0-nrm << endl;
      nrm = 1.0 / nrm;
      dscal(&ndim, &nrm, vec, &one);
    }  
  }
}


///////////////////////////////////////////////////////////////////////////////
///
///  Gram-Schmidt orthonormalize a set of vectors
///
///////////////////////////////////////////////////////////////////////////////
void GramSchmidt(int ndim, int nvec, double *vec, int verbose)
{
  const double RepeatIt = 2;
  int one = 1;
  if (verbose > 5)
    cout << "  Orthonormalize a set of " << nvec << " vectors\n";
  for (int ivec = 0; ivec < nvec; ++ivec) {
    if (verbose > 5)
    cout << "  Vector " << ivec << "\n";
    double *vec_i = vec + ndim*ivec;
    // normalize 
    double nrm = 0;
    nrm = 1.0 / dnrm2(&ndim, vec_i, &one);
    dscal(&ndim, &nrm, vec_i, &one); 
    for (int jvec = 0; jvec < ivec; ++jvec) {
      double *vec_j = vec + ndim*jvec;
      // orthogonalize (twice)
      for (int ri = 0; ri < RepeatIt; ++ri) {
	if (verbose > 5)
	  cout << "  Orthogonalize wrt vector " << jvec << "\n";
	double ovl = -ddot(&ndim, vec_i, &one, vec_j, &one);
	if (verbose > 5)
	  cout << "    overlap " << ri << " = " << fabs(ovl) << "\n";
	daxpy(&ndim, &ovl, vec_j, &one, vec_i, &one);
	nrm = dnrm2(&ndim, vec_i, &one);
	if (verbose > 5)
	  cout << "    norm-loss-" << ri << "  = " << 1-nrm << endl;
	nrm = 1.0 / nrm;
	dscal(&ndim, &nrm, vec_i, &one);
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////////
/// simple Jacobi-Davidson correction vector 
///
///  \f[
///  c = -M r + f M u
///  \f]
///
///  \f$M\f$ is the inverse of \f$\mathrm{diag}(H)-\lambda\f$ as used in standard Davidson
///  \f$r\f$ is the residual vector, \f$u\f$ is the current ritz vector
///  \f$f = (u,Mr)/(u,Mu)\f$
///
///  \f$rcvec = -M r\f$;  \f$rjvec = M u\f$
///////////////////////////////////////////////////////////////////////////////
void DavidsonJacobiCorrectionVector(int ndim, 
				    double lambda, 
				    double *res_vec,
				    double *jd_vec,
				    double *ritz_vec, 
				    double *diagH) 
{
  double shift = 1e-12;
  for (int k = 0; k < ndim; ++k) {
    double mk = diagH[k] - lambda;
    if (fabs(mk) < shift)
      mk = shift;
    mk = 1.0 / mk;
    res_vec[k] *= -mk;
    jd_vec[k] = mk * ritz_vec[k];
  }
  int one = 1;
  double fjd = ddot(&ndim, ritz_vec, &one, res_vec, &one) / ddot(&ndim, ritz_vec, &one, jd_vec, &one);
  daxpy(&ndim, &fjd, jd_vec, &one, res_vec, &one);
}

///////////////////////////////////////////////////////////////////////////////
///
/// standard Davidson correction vector    
///  \f$c = (\lambda_0-\mathrm{diag}(H))^{-1} r\f$
///
///////////////////////////////////////////////////////////////////////////////
void DavidsonCorrectionVector(int ndim, double lambda, double *res_vec, double *diagH) 
{
  double shift = 1e-12;
  for (int k = 0; k < ndim; ++k) {
    double div = lambda - diagH[k];
    if (fabs(div) < shift)
      div = shift;
    res_vec[k] /= div;
  }
}
