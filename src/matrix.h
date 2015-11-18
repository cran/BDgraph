#ifndef matrix_H
#define matrix_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <vector>        // for using vector

using namespace std;

extern "C" {
	void copyMatrix( double A[], double copyA[], int *pxp );

	void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p  );

	void sub_row_mins( double A[], double sub_A[], int *sub, int *p );

	void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p );

	void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p );

	void sub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p );

	void multiply_matrix( double A[], double B[], double C[], int *p_i, int *p_j, int *p_k );

	void inverse( double A[], double A_inv[], int *p );

	void inverse_2x2( double B[], double B_inv[] );

	void cholesky( double A[], double U[], int *p );
	
	void select_edge( long double rates[], int *index_selected_edge, long double *sum_rates, int *qp );
	
	void select_multi_edges( long double rates[], int index_selected_edges[], int *size_index, long double *sum_rates, int *multi_update, int *qp );
}

#endif
