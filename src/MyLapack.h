#ifndef R_MYLAPACK_H
#define R_MYLAPACK_H

#include "util.h"

#ifdef	__cplusplus
extern "C" {
#endif

// Never defined by R itself.
#ifndef La_extern
#define La_extern extern
#endif

// Function for solving Hermitian matrix
int zpotrs( char *uplo, int *n, int *nrhs, Rcomplex *a, int *lda, Rcomplex *b, int *ldb, int *info );

// The Cholesky decomposition for Hermitian matrix
La_extern void
F77_NAME(zpotrf)( const char* uplo, const int* n, Rcomplex* a, const int* lda, int* info );

#ifdef	__cplusplus
}
#endif

#endif /* R_MYLAPACK_H */

