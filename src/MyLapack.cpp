#include "MyLapack.h"

int zpotrs( char *uplo, int *n, int *nrhs, Rcomplex *a, int *lda, Rcomplex *b, int *ldb, int *info )
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   

    Purpose   
    =======   
    ZPOTRS solves a system of linear equations A*X = B with a Hermitian   
    positive definite matrix A using the Cholesky factorization   
    A = U**H*U or A = L*L**H computed by ZPOTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) COMPLEX*16 array, dimension (LDA,N)   
            The triangular factor U or L from the Cholesky factorization   
            A = U**H*U or A = L*L**H, as computed by ZPOTRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
    =====================================================================   
       Test the input parameters.   

       Parameter adjustments */
    /* Table of constant values */
    static Rcomplex c_b1 = {1.,0.};
    
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset;
    /* Local variables */
    static int upper;

    a_dim1   = *lda;
    a_offset = 1 + a_dim1 * 1;
    a       -= a_offset;
    b_dim1   = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b       -= b_offset;

    /* Function Body */
    *info    = 0;
    char up  = 'U';
    upper    = F77_NAME(lsame)( uplo, &up );
    char low = 'L';
    
    if( ! upper && ! F77_NAME(lsame)( uplo, &low ) ) 
    {
		*info = -1;
    }else if( *n < 0 ){
		*info = -2;
    }else if( *nrhs < 0 ){
		*info = -3;
    }else if( *lda < std::max( 1, *n ) ){
		*info = -5;
    }else if( *ldb < std::max( 1, *n ) ){
		*info = -7;
    }
    
    if( *info != 0 )            return 0;
    if( *n == 0 || *nrhs == 0 ) return 0;

	char left_m = 'L', ct = 'C', nt = 'N', nu = 'N';
	
    if( upper ) 
    {
		// Solve A*X = B where A = U'*U.   
		// Solve U'*X = B, overwriting B with X. 
		F77_NAME(ztrsm)( &left_m, &up, &ct, &nu, n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb );

		// Solve U*X = B, overwriting B with X.
		F77_NAME(ztrsm)( &left_m, &up, &nt, &nu, n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb );
    }else{
		// Solve A*X = B where A = L*L'.   
		// Solve L*X = B, overwriting B with X.
		F77_NAME(ztrsm)( &left_m, &low, &nt, &nu, n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb );

		// Solve L'*X = B, overwriting B with X.
		F77_NAME(ztrsm)( &left_m, &low, &ct, &nu, n, nrhs, &c_b1, &a[a_offset], lda, &b[b_offset], ldb );
    }

    return 0;
} /* zpotrs_ */
