#ifndef R_MYLAPACK_H
#define R_MYLAPACK_H

#include <limits>   
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector

#include <R_ext/RS.h>		/* for F77_... */
#include <R_ext/Complex.h>	/* for Rcomplex */
#include <R_ext/BLAS.h>
#include <math.h>

#ifdef	__cplusplus
extern "C" {
#endif

// Never defined by R itself.
#ifndef La_extern
#define La_extern extern
#endif

// Utilities for Lapack-using packages :
// ------------------------------------

/* matrix norms: converting typstr[]  to one of {'M', 'O', 'I', 'F'}
 * or signal error(): */
// La_extern char La_norm_type(const char *typstr);

/* matrix (reciprocal) condition numbers: convert typstr[]  to 'O'(ne) or 'I'(nf)
 * or signal error(): */
// La_extern char La_rcond_type(const char *typstr);


// Function for solving Hermitian matrix
int zpotrs( char *uplo, int *n, int *nrhs, 
	        Rcomplex *a, int *lda, 
	        Rcomplex *b, int *ldb, int *info );

// The Cholesky decomposition for Hermitian matrix
La_extern void
F77_NAME(zpotrf)(const char* uplo, const int* n,
		 Rcomplex* a, const int* lda, int* info);

#ifdef	__cplusplus
}
#endif

#endif /* R_MYLAPACK_H */

// Local variables: ***
// mode: outline-minor ***
// outline-regexp: "^\^L\\|^//[*]+" ***
// End: ***
