#include <boost/fortran/prototype.hpp>

#define FUN BOOST_FORTRAN_FUNCTION
#define SUB BOOST_FORTRAN_SUBROUTINE

SUB( dscal, dscal, DSCAL, (int *)(double *)(double *)(int *) );
SUB( dcopy, dcopy, DCOPY, (int *)(double *)(int *)(double *)(int *)); 
SUB( daxpy, daxpy, DAXPY, (int *)(double *)(double *)(int *)(double *)(int *));

FUN( double, dnrm2, dnrm2, DNRM2, (int *)(double *)(int *) );
FUN( double, ddot, ddot, DDOT, (int *)(const double *)(int *)(const double *)(int *));
SUB( dsyev, dsyev, DSYEV, (const char *)(const char *)(int *)(double *)(int *)(double *)(double *)(int *)(int *));
SUB( dsteqr, dsteqr, DSTEQR, (const char *)(int *)(double *)(double *)(double *)(int *)(double *)(int *));

SUB( dgemv, dgemv, DGEMV, 
    (const char *)(const int *)(const int *)(const double *)(const double *)(const int *)(const double *)
    (const int *)(const double *)(double *)(const int *) );

SUB( dsymv, dsymv, DSYMV, 
     (const char*) // uplo
     (const int&)  // n
     (const double&) // alpha
     (const double *) // ap
     (const int&) // lda
     (const double*) // x
     (const int&) //incx
     (const double&) // beta
     (double *) // y
     (const int&) // incy
     );


SUB( dgemm, dgemm, DGEMM,   // C := alpha * A*B + beta * C
     (const char*) // transA  N or T
     (const char*) // transB  N or T
     (const int&)  // M: rows of C
     (const int&)  // N: columns of C
     (const int&)  // K: columns of op(A) and rows of op(B) 
     (const double&) // alpha
     (const double *) // A
     (const int&) // lda
     (const double*) // B
     (const int&) // ldb
     (const double&) // beta
     (double*) // C
     (const int&) // ldc
     );


SUB( dsytrf, dsytrf, DSYTRF,
     (const char*) // uplo, 
     (int*)        // n, 
     (double*) // a, 
     (int*) // lda, 
     (int*) // ipiv, 
     (double*) // w,
     (int*) // lwork 
     (int*) // info);
     );

SUB ( dsytri, dsytri, DSYTRI,
      (const char*) // uplo, 
      (int*) // n, 
      (double*) // a, 
      (int*) // lda, 
      (int*) // ipiv, 
      (double*) // w, 
      (int*) // info);
      );

SUB( dspmv, dspmv, DSPMV, 
    (const char*) // uplo
    (const int&)  // n
    (const double&) // alpha
    (const double*) // ap
    (const double*) // x
    (const int&) //incx
    (const double&) // beta
    (double *) // y
    (const int&) // incy
    );

//asymmetric matrices 
SUB( dgetrf, dgetrf, DGETRF,
     (int*)        // n,
     (int*)        // n,
     (double*) // a,
     (int*) // lda,
     (int*) // ipiv,
     (int*) // info);
     );

SUB( dgetri, dgetri, DGETRI,
     (int*)        // n,
     (double*) // a,
     (int*) // lda,
     (int*) // ipiv,
     (double*) // w, 
     (int*) // lwork 
     (int*) // info);
     );




#undef SUB
#undef FUN
