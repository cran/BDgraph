#ifndef matrix_H
#define matrix_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

extern "C" {
	void copyMatrix( double A[], double copyA[], int *pxp );

	void subMatrix( double A[], double subA[], int sub[], int *p_sub, int *p  );

	void subRowMins( double A[], double subA[], int *sub, int *p );

	void subRowsMins( double A[], double subA[], int *row, int *col, int *p );

	void subMatrices1( double A[], double A12[], double A22[], int *sub, int *p );

	void subMatrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p );

	void multiplyMatrix( double A[], double B[], double C[], int *p_i, int *p_j, int *p_k );

	void inverse( double A[], double A_inv[], int *p );

	void inverse2x2( double B[], double B_inv[] );

	void cholesky( double A[], double U[], int *p );
}

#endif
