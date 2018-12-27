int DavidsonWorkSize(int ndim, int maxsub, int nroots, int jdflag);
int Davidson(int ndim,        // dimension of H, length of vectors in B and Z
	     int maxsub,      // max subspace dimension
	     int nroots,      // no of roots of H to be found
	     int maxmacro,    // max no of macro iterations
	     int tol,         // norm of ritz vectors must be below 10^-tol
	     int corrflag,    // correction vector flag (see below)
	     int verbose,     // verbocity level
             double *eval,    // returns converged eigenvalues
	     int &nConv,      // returns number of converged eigenpairs
	     double *B,       // space for basis set build by Davidson iteration (contains start vectors on input)
	     double *Z,       // space for H*B vectors
	     double *diag,    // diag(H) for Davidson correction vectors
	     double *davwork, // work space for subspace arrays and ritz vectors 
	     int *inout);      // communication codes with the reverse interface


