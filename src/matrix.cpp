#include "matrix.h"
   
// Takes square matrix A (p x p) and retrieves square sub_matrix B (p_sub x p_sub), dictated by vector sub
void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p  )
{
	int ixp, subixp, psub = *p_sub, pdim = *p;
	
	for( int i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
		
		for( int j = 0; j < psub; j++ )
			sub_A[ixp + j] = A[subixp + sub[j]]; 
	}
}

// Takes square matrix A (p x p) and retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int i, subj = *sub, pdim = *p, subxp = subj * pdim;

	//~ for( i = 0; i < subj; i++ ) sub_A[i] = A[subxp + i];	
	memcpy( sub_A,        A + subxp,            sizeof( double ) * subj );	
	
	//~ for( i = subj + 1; i < pdim; i++ ) sub_A[i - 1] = A[subxp + i];	
	memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// Takes square matrix A (p x p) and retrieves sub_matrix sub_A(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY  FOR SYMMETRIC MATRICES
void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col;
	int sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[l++] = A[sub0p + i]; 
		sub_A[l++] = A[sub1p + i]; 
	}
}

// Takes symmatric matrix A (p x p) and retrieves A_jj, A12(1x(p-1)), A21((p-1)x1), and A22((p-1)x(p-1))
// Like A11=A[j, j], A12=A[j, -j], and A22=A[-j, -j] in R
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
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
void sub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
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
void multiply_matrix( double A[], double B[], double C[], int *p_i, int *p_j, int *p_k )
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
void inverse_2x2( double B[], double B_inv[] )
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
  
// To select an edge for BDMCMC algorithm  
void select_edge( long double rates[], int *index_selected_edge, long double *sum_rates, int *qp )
{
	long double random_value;
	int qp_star = *qp;

	vector<long double> sum_sort_rates( qp_star );
	sum_sort_rates[0] = rates[0];
	for ( int i = 1; i < qp_star; i++ )
		sum_sort_rates[i] = sum_sort_rates[ i - 1 ] + rates[i];
	
	*sum_rates = sum_sort_rates[qp_star - 1];
	random_value = *sum_rates * runif( 0, 1 );

	// To start, find the subscript of the middle position.
	int position;
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	position = upper_bound / 2;      // ( lower_bound + upper_bound ) / 2;

	while( upper_bound - lower_bound > 1 )
	{
		if ( sum_sort_rates[position] > random_value )    
			upper_bound = position;    
		else                                                
			lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	//~ if ( sum_sort_rates[position] < random_value ) position = position + 1;
	*index_selected_edge = ( sum_sort_rates[position] < random_value ) ? ++position : position;
} 
    
// To simultaneously select multiple edges for BDMCMC algorithm  
void select_multi_edges( long double rates[], int index_selected_edges[], int *size_index, long double *sum_rates, int *multi_update, int *qp )
{
	int qp_star = *qp, lower_bound, upper_bound, position;
	long double max_bound, random_value;

	vector<long double> sum_sort_rates( qp_star );
	sum_sort_rates[0] = rates[0];
	for ( int i = 1; i < qp_star; i++ )
		sum_sort_rates[i] = sum_sort_rates[ i - 1 ] + rates[i];
	
	max_bound = sum_sort_rates[qp_star - 1];
	
// ---------- for first edge ---------------------------------------
	// To start, find the subscript of the middle position.
	lower_bound = 0;
	upper_bound = qp_star - 1;
	position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

	random_value = max_bound * runif( 0, 1 );

	while( upper_bound - lower_bound > 1 )
	{
		if ( sum_sort_rates[position] > random_value )    
			upper_bound = position;    
		else                                                
			lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	if ( sum_sort_rates[position] < random_value ) ++position;
	index_selected_edges[0] = position;
// ---------------------------------------------------------------------

	int counter = 1;
	int same;
	for ( int it = 0; it < 200 * *multi_update; it++ )
	{
		if ( counter == *multi_update ) break;
		
		random_value = max_bound * runif( 0, 1 );
	
		// To start, find the subscript of the middle position.
		lower_bound = 0;
		upper_bound = qp_star - 1;
		position   = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

		while( upper_bound - lower_bound > 1 )
		{
			if ( sum_sort_rates[position] > random_value )    
				upper_bound = position;    
			else                                                
				lower_bound = position;     
			
			position = ( lower_bound + upper_bound ) / 2;
		}
		
		if ( sum_sort_rates[position] < random_value ) ++position;
		
		same = 0;
		for ( int i = 0; i < counter; i++ )
			if( index_selected_edges[i] == position )
				++same;

		if ( same == 0 ) index_selected_edges[counter++] = position;
	}
	
	*size_index = counter;
	*sum_rates  = max_bound;
} 
      




  
  
    
