#include "matrix.h"
   
// Takes square matrix A (p x p) and retrieves square submatrix B (p_sub x p_sub), dictated by vector sub
void subMatrix( double A[], double subA[], int sub[], int *p_sub, int *p  )
{
	int ixp, subixp;
	
	for( int i = 0, psub = *p_sub, pdim = *p; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
		
		for( int j = 0; j < psub; j++ )
			subA[ixp + j] = A[subixp + sub[j]]; 
	}
}

// Takes square matrix A (p x p) and retrieves vector subA which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
void subRowMins( double A[], double subA[], int *sub, int *p )
{
	int i, dimSub = *sub, pdim = *p;
	int subxp = dimSub * pdim;

	//~ for( i = 0; i < dimSub; i++ ) subA[i] = A[subxp + i];	
	memcpy( subA,          A + subxp,              sizeof( double ) * dimSub );	
	
	//~ for( i = dimSub + 1; i < pdim; i++ ) subA[i - 1] = A[subxp + i];	
	memcpy( subA + dimSub, A + subxp + dimSub + 1, sizeof( double ) * ( pdim - dimSub - 1 ) );	
}
   
// Takes square matrix A (p x p) and retrieves submatrix subA(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY  FOR SYMMETRIC MATRICES
void subRowsMins( double A[], double subA[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col;
	int sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		subA[l++] = A[sub0p + i]; 
		subA[l++] = A[sub1p + i]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		subA[l++] = A[sub0p + i]; 
		subA[l++] = A[sub1p + i]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		subA[l++] = A[sub0p + i]; 
		subA[l++] = A[sub1p + i]; 
	}
}

// Takes symmatric matrix A (p x p) and retrieves A_jj, A12(1x(p-1)), A21((p-1)x1), and A22((p-1)x(p-1))
// Like A11=A[j, j], A12=A[j, -j], and A22=A[-j, -j] in R
void subMatrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, j, ixpdim, ij, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim;

	memcpy( A12,        A + subxp,            sizeof( double ) * psub );	
	memcpy( A12 + psub, A + subxp + psub + 1, sizeof( double ) * ( pdim - psub - 1 ) );	

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		
		for( j = 0; j < psub; j++ )
			A22[j * p1 + i] = A[ixpdim + j];

		for( j = psub + 1; j < pdim; j++ )
		{
			ij = ixpdim + j;
			A22[(j - 1) * p1 + i] = A[ij];
			A22[i * p1 + j - 1]   = A[ij];
		}
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
			
		for( j = psub + 1; j < pdim; j++ )
			A22[(j - 1) * p1 + i - 1] = A[ixpdim + j];
	}
}

// Takes square matrix A (p x p) and retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
// Like A11=A[e, e], A12=A[e, -e], and A22=A[-e, -e] in R
void subMatrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, j, ixp, ij, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[0] = A[sub0 * pdim + sub0];
	A11[1] = A[sub0 * pdim + sub1];
	A11[2] = A11[1];                   // for symmetric matrices
	A11[3] = A[sub1 * pdim + sub1];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[i + i]     = A[ixp + sub0];
		A12[i + i + 1] = A[ixp + sub1];
	
		for( j = 0; j < sub0; j++ )
			A22[j * p2 + i] = A[ixp + j];

		for( j = sub0 + 1; j < sub1; j++ )
		{
			ij = ixp + j;
			A22[(j - 1) * p2 + i] = A[ij];
			A22[i * p2 + j - 1]   = A[ij];
		}
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[(j - 2) * p2 + i] = A[ij];
			A22[i * p2 + j - 2]   = A[ij];
		}
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		
		A12[i + i - 2] = A[ixp + sub0];
		A12[i + i - 1] = A[ixp + sub1];
	
		for( j = sub0 + 1; j < sub1; j++ )
			A22[(j - 1) * p2 + i - 1] = A[ixp + j];
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[(j - 2) * p2 + i - 1] = A[ij];
			A22[(i - 1) * p2 + j - 2] = A[ij];
		}
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
				
		A12[i + i - 4]     = A[ixp + sub0];
	  //A12[(i - 2) * 2 + 1] = A[i * pdim + sub1];
		A12[i + i - 3]       = A[ixp + sub1];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[(j - 2) * p2 + i - 2] = A[ixp + j];
	}
}
   
////////////////////////////////////////////////////////////////////////////////
//  Multiplies (p_i x p_k) matrix by (p_k x p_j) matrix to give (p_i x p_j) matrix
//  C := A %*% B
void multiplyMatrix( double A[], double B[], double C[], int *p_i, int *p_j, int *p_k )
{
	double alpha = 1.0, beta  = 0.0;
	char trans   = 'N';																	
	F77_NAME(dgemm)( &trans, &trans, p_i, p_j, p_k, &alpha, A, p_i, B, p_k, &beta, C, p_i );
}

// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
void inverse( double A[], double A_inv[], int *p )
{
	int info, dim = *p;
	char uplo = 'U';

	// creating an identity matrix

	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < dim; j++ )
			A_inv[j * dim + i] = (i == j);
	
	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info );
}

// inverse function for symmetric (2 x 2)
void inverse2x2( double B[], double B_inv[] )
{
	double detB = B[0] * B[3] - B[1] * B[1];
	B_inv[0]    = B[3] / detB;
	B_inv[1]    = - B[1] / detB;
	B_inv[2]    = B_inv[1];
	B_inv[3]    = B[0] / detB;
}
    
// Cholesky decomposition of symmetric positive-definite matrix
// WARNING: Matrix you pass is overwritten with the result
// A = U' %*% U
void cholesky( double A[], double U[], int *p )
{
	char uplo = 'U';
	int j, info, dim = *p, i, pxp = dim * dim;
	
	//~ for( i = 0; i < pxp; i++ ) U[i] = A[i]; 
	memcpy( U, A, sizeof( double ) * pxp );	
	
	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info );	
	
	for( i = 0; i < dim; i++ )
		for( j = 0; j < i; j++ )
			U[j * dim + i] = 0.0;
}
    
