// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//    Copyright (C) 2012 - 2022  Reza Mohammadi                                |
//                                                                             |
//    This file is part of BDgraph package.                                    |
//                                                                             |
//   BDgraph is a free software: you can redistribute it and/or modify it      |
//   under the terms of the GNU General Public License as published by the Free|
//   Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>. |
//                                                                             |
//   Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                           |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  
#include "matrix.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves square sub_matrix B (p_sub x p_sub), dictated by vector sub
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
	
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[ i ] * pdim;
		
		for( j = 0; j < psub; j++ )
			sub_A[ ixp + j ] = A[ subixp + sub[ j ] ]; 
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmetric matrix A (p x p) and 
// retrieves upper part of sub_matrix B (p_sub x p_sub), dictated by vector sub
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrix_upper( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
			
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[ i ] * pdim;
			
		for( j = 0; j <= i; j++ )
			sub_A[ ixp + j ] = A[ subixp + sub[ j ] ]; 
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[ j, -j ] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int subj = *sub, pdim = *p, subxp = subj * pdim;

	memcpy( sub_A       , A + subxp           , sizeof( double ) * subj );		
	memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY FOR SYMMETRIC MATRICES
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[ l++ ] = A[ sub0p + i ]; 
		sub_A[ l++ ] = A[ sub1p + i ]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[ l++ ] = A[ sub0p + i ]; 
		sub_A[ l++ ] = A[ sub1p + i ]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[ l++ ] = A[ sub0p + i ]; 
		sub_A[ l++ ] = A[ sub1p + i ]; 
	}
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(p-2 x 2) which is sub cols of matrix A, minus two elements
// Likes A[-(i,j), (i,j)] in R 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_cols_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int subi = *row, subj = *col, pdim = *p, p2 = pdim - 2, subixp = subi * pdim, subjxp = subj * pdim;

	memcpy( sub_A           , A + subixp           , sizeof( double ) * subi );		
	memcpy( sub_A + subi    , A + subixp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );	
	memcpy( sub_A + subj - 1, A + subixp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	

	memcpy( sub_A + p2           , A + subjxp           , sizeof( double ) * subi );		
	memcpy( sub_A + p2 + subi    , A + subjxp + subi + 1, sizeof( double ) * ( subj - subi - 1 ) );	
	memcpy( sub_A + p2 + subj - 1, A + subjxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes symmatric matrix A (p x p) and 
// retrieves A12(1x(p-1)) and A22((p-1)x(p-1))
// Like A12=A[j, -j], and A22=A[-j, -j] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, ixp1, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;
    int size_psub  = sizeof( double ) * psub;
    int size_mpsub = sizeof( double ) * mpsub;

	memcpy( A12,        A + subxp,            size_psub );	
	memcpy( A12 + psub, A + subxp + psub + 1, size_mpsub );	

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		ixp1   = i * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		ixp1   = ( i - 1 ) * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , size_psub );
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, size_mpsub );
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
// Like A11=A[e, e], A12=A[e, -e], and A22=A[-e, -e] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, j, ixp, ij, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[ 0 ] = A[ sub0 * pdim + sub0 ];
	A11[ 1 ] = A[ sub0 * pdim + sub1 ];
	A11[ 2 ] = A11[ 1 ];                   // for symmetric matrices
	A11[ 3 ] = A[ sub1 * pdim + sub1 ];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[ i + i     ] = A[ ixp + sub0 ];
		A12[ i + i + 1 ] = A[ ixp + sub1 ];
	
		for( j = 0; j < sub0; j++ )
			A22[ j * p2 + i ] = A[ ixp + j ];

		for( j = sub0 + 1; j < sub1; j++ )
		{
			ij = ixp + j;
			A22[ ( j - 1 ) * p2 + i ] = A[ ij ];
			A22[ i * p2 + j - 1     ] = A[ ij ];
		}
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[ ( j - 2 ) * p2 + i ] = A[ ij ];
			A22[ i * p2 + j - 2     ] = A[ ij ];
		}
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		
		A12[ i + i - 2 ] = A[ ixp + sub0 ];
		A12[ i + i - 1 ] = A[ ixp + sub1 ];
	
		for( j = sub0 + 1; j < sub1; j++ )
			A22[ ( j - 1 ) * p2 + i - 1 ] = A[ ixp + j ];
		
		for( j = sub1 + 1; j < pdim; j++ )
		{
			ij = ixp + j;
			A22[ ( j - 2 ) * p2 + i - 1 ] = A[ ij ];
			A22[ ( i - 1 ) * p2 + j - 2 ] = A[ ij ];
		}
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
				
		A12[ i + i - 4 ] = A[ ixp + sub0 ];
		A12[ i + i - 3 ] = A[ ixp + sub1 ];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[ ( j - 2 ) * p2 + i - 2 ] = A[ ixp + j ];
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Takes square matrix A (p x p) and 
// retrieves A11_inv ( 2 x 2 ), A21 ( ( p - 2 ) x 2 ), and A22 ( ( p - 2 ) x ( p - 2 ) )
// Like A11_inv=inv ( A[ e, e ] ), A21 = A[ -e, e ], and A22 = A[ -e, -e ] in R
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int *row, int *col, int *p )
{
	int i, ixp, ixp2, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;
	int sub0xp = sub0 * pdim, sub1xp = sub1 * pdim, sub0_plus = sub0 + 1, sub1_plus = sub1 + 1;
	
	double a11 = A[ sub0 * pdim + sub0 ];
	double a12 = A[ sub0 * pdim + sub1 ];
	double a22 = A[ sub1 * pdim + sub1 ];

	double det_A11 = a11 * a22 - a12 * a12;
	A11_inv[ 0 ]   = a22 / det_A11;
	A11_inv[ 1 ]   = - a12 / det_A11;
	A11_inv[ 2 ]   = A11_inv[ 1 ];
	A11_inv[ 3 ]   = a11 / det_A11;
	
	int size_sub0      = sizeof( double ) * sub0;
	int size_sub1_sub0 = sizeof( double ) * ( sub1 - sub0_plus );
	int size_pdim_sub0 = sizeof( double ) * ( pdim - sub1_plus );
	
	memcpy( A21           , A + sub0xp            , size_sub0 );		
	memcpy( A21 + sub0    , A + sub0xp + sub0_plus, size_sub1_sub0 );	
	memcpy( A21 + sub1 - 1, A + sub0xp + sub1_plus, size_pdim_sub0 );	

	memcpy( A21 + p2           , A + sub1xp            , size_sub0 );		
	memcpy( A21 + p2 + sub0    , A + sub1xp + sub0_plus, size_sub1_sub0 );	
	memcpy( A21 + p2 + sub1 - 1, A + sub1xp + sub1_plus, size_pdim_sub0 );	
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp  = i * pdim;
		ixp2 = i * p2;

		memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );	
	}
 
	for( i = sub0_plus; i < sub1; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 1 ) * p2;

		memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );	
	}
	
	for( i = sub1_plus; i < pdim; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 2 ) * p2;
				
		memcpy( A22 + ixp2           , A + ixp            , size_sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, size_sub1_sub0 );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, size_pdim_sub0 );		
	}
}
      
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void inverse( double A[], double A_inv[], int *p )
{
	int info, dim = *p;
	char uplo = 'U';

	// creating an identity matrix
	#pragma omp parallel for
	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < dim; j++ )
			A_inv[ j * dim + i ] = ( i == j );
	
	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info FCONE );
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for symmetric (2 x 2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void inverse_2x2( double B[], double B_inv[] )
{
	double detB = B[ 0 ] * B[ 3 ] - B[ 1 ] * B[ 1 ];
	B_inv[ 0 ]  = B[ 3 ] / detB;
	B_inv[ 1 ]  = - B[ 1 ] / detB;
	B_inv[ 2 ]  = B_inv[ 1 ];
	B_inv[ 3 ]  = B[ 0 ] / detB;
} 
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Cholesky decomposition of symmetric positive-definite matrix
// A = U' %*% U
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void cholesky( double A[], double U[], int *p )
{
	char uplo = 'U';
	int info, dim = *p;
	
	memcpy( U, A, sizeof( double ) * dim * dim );	
	
	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info FCONE );	
	
	#pragma omp parallel for
	for( int i = 0; i < dim; i++ )
		for( int j = 0; j < i; j++ )
			U[ j * dim + i ] = 0.0;
}
  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//  Determinant of a symmetric positive-definite matrix ( A )
//  > > > > > > > > >  WARNING: Matrix you pass is overwritten < < < < < < < < < 
//  For any symmetric PD Matrix A: |A| = |T| ^ 2, where T is cholesky 
// decomposition of A. Thus, |T| = \prod_{i = 1}^p T_{ii}.
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void determinant( double A[], double *det_A, int *p )
{
	char uplo = 'U';
	int info, dim = *p, dim1 = dim + 1;
	
	F77_NAME(dpotrf)( &uplo, &dim, &A[0], &dim, &info FCONE );	

	double result = 1;
	for( int i = 0; i < dim; i++ ) result *= A[ i * dim1 ];
	
	*det_A = result * result;
}
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To select an edge for BDMCMC algorithm  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_edge( double rates[], int *index_selected_edge, double *sum_rates, int *qp )
{
	int qp_star = *qp;

	// rates = sum_sort_rates
	vector<double>cumulative_rates( qp_star, 0.0 );
	cumulative_rates[ 0 ] = rates[ 0 ];
	for( int i = 1; i < qp_star; i++ )
		cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
	
	*sum_rates = cumulative_rates[ qp_star - 1 ];
	
	// GetRNGstate();
	double random_value = *sum_rates * unif_rand(); // Rf_runif( 0.0, *sum_rates );
	// PutRNGstate();

	//int counter = 0;
	//while( random_value > cumulative_rates[ counter ] )	++counter;
	//*index_selected_edge = counter;
	 
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	int position    = upper_bound / 2;  // ( lower_bound + upper_bound ) / 2;

	while( upper_bound - lower_bound > 1 )
	{
		 //if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	*index_selected_edge = ( cumulative_rates[ position ] < random_value ) ? ++position : position;
} 

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// select_edge for bd_for_ts
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_edge_ts( long double rates[], int *index_selected_edge, long double *sum_rates, int *qp )
{
	int qp_star = *qp;

	// rates = sum_sort_rates
	vector<long double>cumulative_rates( qp_star, 0.0 );
	cumulative_rates[ 0 ] = rates[ 0 ];
	for( int i = 1; i < qp_star; i++ )
		cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
	
	*sum_rates = cumulative_rates[ qp_star - 1 ];
	
	// GetRNGstate();
	long double random_value = *sum_rates * unif_rand();  // Rf_runif( 0.0, *sum_rates );
	// PutRNGstate();

	//int counter = 0;
	//while( random_value > cumulative_rates[ counter ] )	++counter;
	//*index_selected_edge = counter;
	 
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	int position    = upper_bound / 2;  // ( lower_bound + upper_bound ) / 2;

	while( upper_bound - lower_bound > 1 )
	{
		 //if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	*index_selected_edge = ( cumulative_rates[ position ] < random_value ) ? ++position : position;
} 
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To simultaneously select multiple edges for BDMCMC algorithm  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void select_multi_edges( double rates[], int index_selected_edges[], int *size_index, double *sum_rates, int *multi_update, int *qp )
{
	int i, qp_star = *qp, qp_star_1 = qp_star - 1;

	// rates = sum_sort_rates
	vector<double>cumulative_rates( qp_star, 0.0 );
	cumulative_rates[ 0 ] = rates[ 0 ];
	for ( int i = 1; i < qp_star; i++ )
		cumulative_rates[ i ] = cumulative_rates[ i - 1 ] + rates[ i ];
	
	double max_bound = cumulative_rates[ qp_star_1 ];
	
	// - - - - - - for first edge - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star_1;
	int position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

	//GetRNGstate();
	double random_value = max_bound * unif_rand();
	//PutRNGstate();

	while( upper_bound - lower_bound > 1 )
	{
		//if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	if ( cumulative_rates[position] < random_value ) ++position;
	index_selected_edges[0] = position;
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

	int counter = 1, same;
	//GetRNGstate();
	for( int it = 0; it < 200 * *multi_update; it++ )
	{
		if( counter == *multi_update ) break;
		
		random_value = max_bound * unif_rand();
	
		// To start, find the subscript of the middle position.
		lower_bound = 0;
		upper_bound = qp_star_1;
		position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

		while( upper_bound - lower_bound > 1 )
		{
			// if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
			( cumulative_rates[ position ] > random_value ) ? upper_bound = position : lower_bound = position;     
			
			position = ( lower_bound + upper_bound ) / 2;
		}
		
		if( cumulative_rates[position] < random_value ) ++position;
		
		same = 0;
		for( i = 0; i < counter; i++ )
			if( index_selected_edges[ i ] == position )
				++same;

		if( same == 0 ) index_selected_edges[ counter++ ] = position;
	}
	//PutRNGstate();

	*size_index = counter;
	*sum_rates  = max_bound;
} 
         
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Parallel Computation for birth-death rates for BD-MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rates_bdmcmc_parallel( double rates[], double log_ratio_g_prior[], int G[], int index_row[], int index_col[], int *sub_qp, double Ds[], double Dsijj[],
				            double sigma[], double K[], int *b, int *p )
{
	int b1 = *b, one = 1, two = 2, dim = *p, p1 = dim - 1, p2 = dim - 2, dim1 = dim + 1, p2x2 = ( dim - 2 ) * 2;
	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	

	#pragma omp parallel
	{
		int i, j, k, ij, jj, nu_star;
		double Dsjj, sum_diag, K022, a11, sigmajj_inv, log_rate;

		double *K121         = new double[ 4 ];  
		double *Kj12         = new double[ p1 ];  
		double *sigmaj12     = new double[ p1 ];  
		double *sigmaj22     = new double[ p1 * p1 ];  
		double *Kj12xK22_inv = new double[ p1 ];  
		
		double *K21                 = new double[ p2x2 ];  
		double *sigma21             = new double[ p2x2 ];  
		double *sigma22             = new double[ p2 * p2 ];  
		double *sigma11_inv         = new double[ 4 ];  
		double *sigma21xsigma11_inv = new double[ p2x2 ];  
		double *K12xK22_inv         = new double[ p2x2 ];  

		#pragma omp for
		for( int counter = 0; counter < *sub_qp; counter++ )
		{
			i  = index_row[ counter ];
			j  = index_col[ counter ];
			ij = j * dim + i;
			jj = j * dim1;
			
			Dsjj = Ds[ jj ];
			
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
			sigmajj_inv = - 1.0 / sigma[ jj ];
			F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 FCONE );
			
			// For (i,j) = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
			sub_row_mins( &K[0], &Kj12[0], &j, &dim );  // Kj12 = K[j, -j]  
			Kj12[ i ] = 0.0;                            // Kj12[1,i] = 0

			// Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
			F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one FCONE );
			
			// K022 = Kj12xK22_inv %*% t(Kj12)
			K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			

			// For (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
			sub_cols_mins( &K[0], &K21[0], &i, &j, &dim );  // K21 = K[-e, e]  
			
			sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], &i, &j, &dim );

			// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
			F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 FCONE FCONE );

			// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
			F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 FCONE FCONE );

			// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
			F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two FCONE FCONE );  
			
			// K121 = K12xK22_inv %*% K21													
			F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two FCONE FCONE );		
			// Finished (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

			a11      = K[ i * dim1 ] - K121[ 0 ];	
			sum_diag = Dsjj * ( K022 - K121[ 3 ] ) - Ds[ ij ] * ( K121[ 1 ] + K121[ 2 ] );

			// nu_star = b + sum( Gf[,i] * Gf[,j] )
			nu_star = b1;
			//for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k];   
			for( k = 0; k < dim; k++ ) nu_star += G[ i * dim + k ] * G[ j * dim + k ];   
			//nu_star = F77_NAME(ddot)( &dim, &G[0] + ixdim, &one, &G[0] + jxdim, &one );
			nu_star = 0.5 * nu_star;

			log_rate = ( G[ ij ] )   
				? 0.5 * log( 2.0 * Dsjj / a11 ) + lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsijj[ ij ] * a11 + sum_diag )
				: 0.5 * log( 0.5 * a11 / Dsjj ) - lgammafn( nu_star + 0.5 ) + lgammafn( nu_star ) + 0.5 * ( Dsijj[ ij ] * a11 + sum_diag );
			
			//log_rate = ( G[ij] ) ? log_rate - log( static_cast<double>( g_prior[ij] / ( 1 - g_prior[ij] ) ) ) : log_rate + log( static_cast<double>( g_prior[ij] / ( 1 - g_prior[ij] ) ) );
			log_rate = ( G[ ij ] ) ? log_rate - log_ratio_g_prior[ ij ] : log_rate + log_ratio_g_prior[ ij ];

			rates[ counter ] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
		}
		delete[] K121;  
		delete[] Kj12;  
		delete[] sigmaj12;  
		delete[] sigmaj22;  
		delete[] Kj12xK22_inv;  
		
		delete[] K21;  
		delete[] sigma21;  
		delete[] sigma22;  
		delete[] sigma11_inv;  
		delete[] sigma21xsigma11_inv;  
		delete[] K12xK22_inv;  
	}
}
     
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Parallel Computation for birth-death rates for complex BD-MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rates_cbdmcmc_parallel( long double log_rates[], double log_ratio_g_prior[], int G[], int index_row[], int index_col[], int *sub_qp,
				             double r_Ds[], double i_Ds[], double r_sigma[], double i_sigma[], double r_K[], double i_K[], int *b, int *p )
{
	int b1 = *b, one = 1, two = 2, dim = *p, p1 = dim - 1, p2 = dim - 2, dim1 = dim + 1, p2x2 = p2 * 2, p2xp2 = p2 * p2;
	double alpha = 1.0, beta = 0.0, dmone = -1.0;
	char transT = 'T', transN = 'N';																	

	#pragma omp parallel
	{
		int i, j, k, ij, jj, rowCol, nu_star;
		double r_Dsjj, r_Dsij, i_Dsij, sum_diag, r_K022, i_K022, a11, r_sigmaj11, i_sigmaj11;
		double coef, epower, I_const, log_rate;

		double *r_K121     = new double[ 4 ];  
		double *i_K121     = new double[ 4 ];  
		double *r_Kj12     = new double[ p1 ];        //K[j,-j]
		double *i_Kj12     = new double[ p1 ];
		double *r_sigmaj12 = new double[ p1 ];        //sigma[j,-j]
		double *i_sigmaj12 = new double[ p1 ];
		double *r_sigmaj22 = new double[ p1 * p1 ];   //sigma[-j,-j]
		double *i_sigmaj22 = new double[ p1 * p1 ];
		double *r_Kj22_inv = new double[ p1 * p1 ]; 
		double *i_Kj22_inv = new double[ p1 * p1 ]; 
		double *i12xi22_j  = new double[ p1 ];
		double *r12xi22_j  = new double[ p1 ];
		double *i12xi22    = new double[ p2x2 ];
		double *r12xi22    = new double[ p2x2 ];
		double *r21xr11    = new double[ p2x2 ];
		double *i21xr11    = new double[ p2x2 ];

		double *r_K12         = new double[ p2x2 ];  //K[e,-e]
		double *i_K12         = new double[ p2x2 ]; 
		double *r_sigma11     = new double[ 4 ];     //sigma[e,e]
		double *i_sigma11     = new double[ 4 ];
		double *r_sigma12     = new double[ p2x2 ];  //sigma[e,-e]
		double *i_sigma12     = new double[ p2x2 ]; 
		double *r_sigma22     = new double[ p2xp2 ]; //sigma[-e,-e]
		double *i_sigma22     = new double[ p2xp2 ];  
		double *r_sigma11_inv = new double[ 4 ]; 
		double *i_sigma11_inv = new double[ 4 ]; 
		double *r_sigma2112   = new double[ p2xp2 ];
		double *i_sigma2112   = new double[ p2xp2 ];  
		double *r_K22_inv     = new double[ p2xp2 ];
		double *i_K22_inv     = new double[ p2xp2 ];
		double *r_K12xK22_inv = new double[ p2x2 ];  
		double *i_K12xK22_inv = new double[ p2x2 ];  

		#pragma omp for
		for( int counter = 0; counter < *sub_qp; counter++ )
		{
			i = index_row[ counter ];
			j = index_col[ counter ];

			jj     = j * dim1;
			r_Dsjj = r_Ds[jj];
			
			r_sigmaj11 = r_sigma[ jj ];        // sigma[j, j]  
			i_sigmaj11 = i_sigma[ jj ]; 			
			sub_matrices1( &r_sigma[0], &r_sigmaj12[0], &r_sigmaj22[0], &j, &dim ); // sigmaj22 = sigma[-j,-j]
			Hsub_matrices1( &i_sigma[0], &i_sigmaj12[0], &i_sigmaj22[0], &j, &dim );

			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj21 * sigmaj12 / sigmaj11
			for( int row = 0; row < p1; row++ )       
				for( int col = 0; col < p1; col++ )
				{
					rowCol               = col * p1 + row;
					r_Kj22_inv[ rowCol ] = r_sigmaj22[rowCol] - (r_sigmaj11*(r_sigmaj12[row]*r_sigmaj12[col] + i_sigmaj12[row]*i_sigmaj12[col]) + i_sigmaj11*(r_sigmaj12[row]*i_sigmaj12[col] - i_sigmaj12[row]*r_sigmaj12[col]))/(r_sigmaj11*r_sigmaj11 + i_sigmaj11*i_sigmaj11);
					i_Kj22_inv[ rowCol ] = i_sigmaj22[rowCol] - (r_sigmaj11*(r_sigmaj12[row]*i_sigmaj12[col] - i_sigmaj12[row]*r_sigmaj12[col]) - i_sigmaj11*(r_sigmaj12[row]*r_sigmaj12[col] + i_sigmaj12[row]*i_sigmaj12[col]))/(r_sigmaj11*r_sigmaj11 + i_sigmaj11*i_sigmaj11);
				}		
			
			ij     = j * dim + i;
			r_Dsij = r_Ds[ ij ];
			i_Dsij = i_Ds[ ij ];

			// For (i,j) = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
			sub_row_mins( &r_K[0], &r_Kj12[0], &j, &dim );   // Kj12 = K[j, -j]
			Hsub_row_mins( &i_K[0], &i_Kj12[0], &j, &dim );  
			r_Kj12[ i ] = 0.0;                             // Kj12[1,i] = 0
			i_Kj12[ i ] = 0.0;

			// Kj12xK22_inv = Kj12 %*% Kj22_inv
			F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &i_Kj12[0], &one, &i_Kj22_inv[0], &p1, &beta, &i12xi22_j[0], &one FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &r_Kj12[0], &one, &r_Kj22_inv[0], &p1, &dmone, &i12xi22_j[0], &one FCONE FCONE );				
			F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &r_Kj12[0], &one, &i_Kj22_inv[0], &p1, &beta, &r12xi22_j[0], &one FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &one, &p1, &p1, &alpha, &i_Kj12[0], &one, &r_Kj22_inv[0], &p1, &alpha, &r12xi22_j[0], &one FCONE FCONE );				
			// K022  <- Kj12 %*% solve( K0[-j, -j] ) %*% t(Kj12) = c
			F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &r12xi22_j[0], &one, &i_Kj12[0], &p1, &beta, &r_K022, &one FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &i12xi22_j[0], &one, &r_Kj12[0], &p1, &alpha, &r_K022, &one FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &i12xi22_j[0], &one, &i_Kj12[0], &p1, &beta, &i_K022, &one FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &one, &one, &p1, &alpha, &r12xi22_j[0], &one, &r_Kj12[0], &p1, &dmone, &i_K022, &one FCONE FCONE );

			// For (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
			sub_rows_mins( &r_K[0], &r_K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
			Hsub_rows_mins( &i_K[0], &i_K12[0], &i, &j, &dim );  // K12 = K[e, -e]  
			
			sub_matrices( &r_sigma[0], &r_sigma11[0], &r_sigma12[0], &r_sigma22[0], &i, &j, &dim ); //r_sigma[e,e], r_sigma[e,-e], r_sigma[-e,-e]
			Hsub_matrices( &i_sigma[0], &i_sigma11[0], &i_sigma12[0], &i_sigma22[0], &i, &j, &dim ); //i_sigma[e,e], i_sigma[e,-e], i_sigma[-e,-e]

			// solve( sigma[e, e] )
			cinverse_2x2( &r_sigma11[0], &i_sigma11[0], &r_sigma11_inv[0], &i_sigma11_inv[0] );
			
			// sigma21 %*% sigma11_inv = t(sigma12) %*% sigma11_inv
			F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &r_sigma12[0], &two, &r_sigma11_inv[0], &two, &beta, &r21xr11[0], &p2 FCONE FCONE );
			F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &i_sigma12[0], &two, &i_sigma11_inv[0], &two, &alpha, &r21xr11[0], &p2 FCONE FCONE );
			F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &i_sigma12[0], &two, &r_sigma11_inv[0], &two, &beta, &i21xr11[0], &p2 FCONE FCONE );
			F77_NAME(dgemm)( &transT, &transN, &p2, &two, &two, &alpha, &r_sigma12[0], &two, &i_sigma11_inv[0], &two, &dmone, &i21xr11[0], &p2 FCONE FCONE );				
			// sigma2112 = sigma21xsigma11_inv %*% sigma12 = sigma[-e,e] %*% solve(sigma[e,e]) %*% sigma[e,-e]
			F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &i21xr11[0], &p2, &i_sigma12[0], &two, &beta, &r_sigma2112[0], &p2 FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &r21xr11[0], &p2, &r_sigma12[0], &two, &dmone, &r_sigma2112[0], &p2 FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &r21xr11[0], &p2, &i_sigma12[0], &two, &beta, &i_sigma2112[0], &p2 FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &p2, &p2, &two, &alpha, &i21xr11[0], &p2, &r_sigma12[0], &two, &alpha, &i_sigma2112[0], &p2 FCONE FCONE );

			// solve( K[-e, -e] ) = sigma22 - sigma2112
			for( k = 0; k < p2xp2 ; k++ ) 
			{
				r_K22_inv[ k ] = r_sigma22[ k ] - r_sigma2112[ k ];
				i_K22_inv[ k ] = i_sigma22[ k ] - i_sigma2112[ k ];
			}

			// K12 %*% K22_inv
			F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &i_K12[0], &two, &i_K22_inv[0], &p2, &beta, &i12xi22[0], &two FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &r_K12[0], &two, &r_K22_inv[0], &p2, &dmone, &i12xi22[0], &two FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &r_K12[0], &two, &i_K22_inv[0], &p2, &beta, &r12xi22[0], &two FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transN, &two, &p2, &p2, &alpha, &i_K12[0], &two, &r_K22_inv[0], &p2, &alpha, &r12xi22[0], &two FCONE FCONE );
			// K121 <- K[e, -e] %*% solve( K[-e, -e] ) %*% t(K[e, -e]) 		
			F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &r12xi22[0], &two, &i_K12[0], &two, &beta, &r_K121[0], &two FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &i12xi22[0], &two, &r_K12[0], &two, &alpha, &r_K121[0], &two FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &i12xi22[0], &two, &i_K12[0], &two, &beta, &i_K121[0], &two FCONE FCONE );
			F77_NAME(dgemm)( &transN, &transT, &two, &two, &p2, &alpha, &r12xi22[0], &two, &r_K12[0], &two, &dmone, &i_K121[0], &two FCONE FCONE );											
			// Finished (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
			
			nu_star = b1;
			for( k = 0; k < dim; k++ ) // nu_star = b + sum( Gf[,i] * Gf[,j] )
				nu_star += G[ i * dim + k ] * G[ j * dim + k ]; 
						
			I_const  = lgammafn( 0.5 * ( nu_star + 1 ) ) - lgammafn( 0.5 * nu_star ); //I

			a11      = r_K[i * dim1] - r_K121[0]; //k_ii - k_ii^1
			sum_diag = r_Dsjj*(r_K022 - r_K121[3]) - (r_Dsij*r_K121[1] - i_Dsij*i_K121[1]) - (r_Dsij*r_K121[2] + i_Dsij*i_K121[2]); //tr(D*(K0-K1))

			coef     = ( r_Dsij * r_Dsij + i_Dsij * i_Dsij ) / r_Dsjj;
			epower   = coef * a11 + sum_diag;			
			log_rate = ( G[ij] ) ? I_const + log( r_Dsjj ) - log( a11 ) - epower : log( a11 ) - log( r_Dsjj ) + epower - I_const;
			
			//log_rates[counter] += log_rate;  // Computer the rate in log space			
			log_rates[counter] += log_rate;
		}
		delete[] r_K121;  
		delete[] i_K121;  
		delete[] r_Kj12;  
		delete[] i_Kj12;  
		delete[] r_sigmaj12;  
		delete[] i_sigmaj12;  
		delete[] r_sigmaj22;  
		delete[] i_sigmaj22;  
		delete[] r_Kj22_inv;
		delete[] i_Kj22_inv;
		delete[] i12xi22_j;
		delete[] r12xi22_j;
		delete[] i12xi22;
		delete[] r12xi22;
		delete[] r21xr11;
		delete[] i21xr11;
		
		delete[] r_K12;  
		delete[] i_K12;  
		delete[] r_sigma11;  
		delete[] i_sigma11;  
		delete[] r_sigma12;  
		delete[] i_sigma12;  
		delete[] r_sigma22;  
		delete[] i_sigma22;  
		delete[] r_sigma11_inv;  
		delete[] i_sigma11_inv;  
		delete[] r_sigma2112;  
		delete[] i_sigma2112;  
		delete[] r_K22_inv;  
		delete[] i_K22_inv;  
		delete[] r_K12xK22_inv;  
		delete[] i_K12xK22_inv;  
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// computing birth/death rate or alpha for element (i,j)
// it is for double Metropolis-Hasting algorithms
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_H_ij( double K[], double sigma[], double *log_Hij, int *selected_edge_i, int *selected_edge_j,
               double Kj12[], double Kj12xK22_inv[], double K12[], double K12xK22_inv[], double K121[], 
               double sigmaj12[], double sigmaj22[], double sigma12[], double sigma22[], double sigma11_inv[], double sigma21xsigma11_inv[],
               int *dim, int *p1, int *p2, int *jj,
               double *Dsijj, double *Dsij, double *Dsjj )
{
	int one = 1, two = 2;
	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	
	
	//double sigmaj11 = sigma[*jj];        // sigma[j, j]  
	sub_matrices1( sigma, sigmaj12, sigmaj22, selected_edge_j, dim );

	// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
	// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
	double sigmajj_inv = - 1.0 / sigma[ *selected_edge_j * ( *dim + 1 ) ];
	F77_NAME(dsyr)( &sideL, p1, &sigmajj_inv, sigmaj12, &one, sigmaj22, p1 FCONE );

	// For (i,j) = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
	sub_row_mins( K, Kj12, selected_edge_j, dim );   // K12 = K[j, -j]  
	Kj12[ *selected_edge_i ] = 0.0;                  // K12[1,i] = 0

	// Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
	F77_NAME(dsymv)( &sideL, p1, &alpha, &sigmaj22[0], p1, Kj12, &one, &beta, Kj12xK22_inv, &one FCONE );
	
	// K022 = Kj12xK22_inv %*% t(Kj12)
	double K022 = F77_NAME(ddot)( p1, Kj12xK22_inv, &one, Kj12, &one );			

	// For (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	sub_cols_mins( K, K12, selected_edge_i, selected_edge_j, dim );   // K21 = K[-e, e] 
	
	sub_matrices_inv( sigma, sigma11_inv, sigma12, sigma22, selected_edge_i, selected_edge_j, dim );

	// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
	F77_NAME(dgemm)( &transN, &transN, p2, &two, &two, &alpha, sigma12, p2, sigma11_inv, &two, &beta, sigma21xsigma11_inv, p2 FCONE FCONE );

	// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
	F77_NAME(dgemm)( &transN, &transT, p2, p2, &two, &alpha1, sigma21xsigma11_inv, p2, sigma12, p2, &beta1, sigma22, p2 FCONE FCONE );

	// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
	F77_NAME(dgemm)( &transT, &transN, &two, p2, p2, &alpha, K12, p2, sigma22, p2, &beta, K12xK22_inv, &two FCONE FCONE );  
	
	// K121 = K12xK22_inv %*% K21													
	F77_NAME(dgemm)( &transN, &transN, &two, &two, p2, &alpha, K12xK22_inv, &two, K12, p2, &beta, K121, &two FCONE FCONE );		
	// Finished (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	double a11      = K[*selected_edge_i * *dim + *selected_edge_i] - K121[0];	
	double sum_diag = *Dsjj * ( K022 - K121[3] ) - *Dsij * ( K121[1] + K121[2] );

	// Dsijj = Dsii - Dsij * Dsij / Dsjj;
	//*log_Hij = ( log( static_cast<double>(*Dsjj) ) - log( static_cast<double>(a11) ) + *Dsijj * a11 - sum_diag ) / 2;
	*log_Hij = 0.5 * ( log( *Dsjj / a11 ) + *Dsijj * a11 - sum_diag );
}    
     
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Parallel Computation for birth-death rates for double BD-MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rates_bdmcmc_dmh_parallel( double rates[], double log_ratio_g_prior[], int G[], int index_row[], int index_col[], int *sub_qp, double Ds[], double D[],
				            double sigma[], double K[], double sigma_dmh[], 
				            double K_dmh[], int *b, int *p )
{
	int dim = *p, p1 = dim - 1, p2 = dim - 2, p2x2 = ( dim - 2 ) * 2;

	#pragma omp parallel
	{
		int index_rate_j, i, j, ij, jj;
		double Dsjj, Dsij, Dsijj, Dij, Dijj, Djj, log_rate;

		double *K121                = new double[ 4 ];  
		double *Kj12                = new double[ p1 ];  
		double *sigmaj12            = new double[ p1 ];  
		double *sigmaj22            = new double[ p1 * p1 ];  
		double *Kj12xK22_inv        = new double[ p1 ];  
		double *K21                 = new double[ p2x2 ];  
		double *sigma12             = new double[ p2x2 ];  
		double *sigma22             = new double[ p2 * p2 ];  
		double *sigma11_inv         = new double[ 4 ];  
		double *sigma21xsigma11_inv = new double[ p2x2 ];  
		double *K12xK22_inv         = new double[ p2x2 ];
		
		double *K12                 = new double[ p2x2 ];
  
		#pragma omp for
		for( j = 1; j < dim; j++ )
		{			
			index_rate_j = ( j * ( j - 1 ) ) / 2;

			jj   = j * dim + j;
			Dsjj = Ds[ jj ];
			Djj  = D[ jj ];

			for( i = 0; i < j; i++ )
			{
				ij    = j * dim + i;
				Dsij  = Ds[ ij ];
				Dsijj = - Dsij * Dsij / Dsjj;
				Dij   = D[ ij ];
				Dijj  = - Dij * Dij / Djj;

				double logH_ij, logI_p;
				
				log_H_ij( &K[0], &sigma[0], &logH_ij, &i, &j,
					   &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
					   &dim, &p1, &p2, &jj,
					   &Dsijj, &Dsij, &Dsjj );

				log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &i, &j,
					   &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
					   &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
					   &dim, &p1, &p2, &jj,
					   &Dijj, &Dij, &Djj );
				
				//log_rate = ( G[ ij ] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );				
				log_rate = ( G[ ij ] ) ? ( logH_ij - logI_p ) - log_ratio_g_prior[ ij ] : ( logI_p - logH_ij ) + log_ratio_g_prior[ ij ];				
				rates[ index_rate_j + i ] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;
			}
		}	
		
		delete[] K121;  
		delete[] Kj12;  
		delete[] sigmaj12;  
		delete[] sigmaj22;  
		delete[] Kj12xK22_inv;  		
		delete[] K21;  
		delete[] sigma12;  
		delete[] sigma22;  
		delete[] sigma11_inv;  
		delete[] sigma21xsigma11_inv;  
		delete[] K12xK22_inv;  
		delete[] K12;  
	}
}
     	
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// NEW for Lang codes for Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int i, l = 0, subj = *sub, pdim = *p, subxp = subj * pdim;

	for( i = 0; i < subj; i++ )
		sub_A[ l++ ] = -A[ subxp + i ];
	
	for( i = subj + 1; i < pdim; i++ )
		sub_A[ l++ ] = -A[ subxp + i ];
}
      
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// For Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[ l++ ] = -A[ sub0p + i ]; 
		sub_A[ l++ ] = -A[ sub1p + i ]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[ l++ ] = -A[ sub0p + i ]; 
		sub_A[ l++ ] = -A[ sub1p + i ]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[ l++ ] = -A[ sub0p + i ]; 
		sub_A[ l++ ] = -A[ sub1p + i ]; 
	}
}
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sub_matrices1 for Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;

	for( i = 0; i < psub; i++ )
		A12[ i ] = -A[ subxp + i ];
	
	for( i = psub; i < pdim - 1; i++ )
		A12[ i ] = -A[ subxp + i + 1 ];

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		memcpy( A22 + i * p1       , A + ixpdim           , sizeof( double ) * psub );
		memcpy( A22 + i * p1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		memcpy( A22 + ( i - 1 ) * p1       , A + ixpdim           , sizeof( double ) * psub);
		memcpy( A22 + ( i - 1 ) * p1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}
}
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sub_matrices for Hermitian matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void Hsub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, i1, i2, ixp, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[ 0 ] = A[ sub0 * pdim + sub0 ];
	A11[ 1 ] = A[ sub0 * pdim + sub1 ];
	A11[ 2 ] = -A11[ 1 ];                   // for symmetric matrices
	A11[ 3 ] = A[ sub1 * pdim + sub1 ];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[ i + i ]     = A[ ixp + sub0 ];
		A12[ i + i + 1 ] = A[ ixp + sub1 ];

		memcpy( A22 + i * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		i1 = i - 1;

		A12[ i + i - 2 ] = A[ ixp + sub0 ];
		A12[ i + i - 1 ] = A[ ixp + sub1 ];

		memcpy( A22 + i1 * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i1 * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i1 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
		i2  = i - 2;
				
		A12[ i + i - 4 ] = A[ ixp + sub0 ];
		A12[ i + i - 3 ] = A[ ixp + sub1 ];

		memcpy( A22 + i2 * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i2 * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i2 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );		
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// inverse function for Hermitian (2 x 2)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void cinverse_2x2( double r_B[], double i_B[], double r_B_inv[], double i_B_inv[] )
{
	double r_det = r_B[0] * r_B[3] - i_B[0] * i_B[3] - ( r_B[1] * r_B[1] + i_B[1] * i_B[1] );
	double i_det = r_B[0] * i_B[3] + i_B[0] * r_B[3];
	double mod   = r_det * r_det + i_det * i_det;
	
	r_B_inv[0] =  ( r_B[3] * r_det + i_B[3] * i_det ) / mod;
	i_B_inv[0] =  ( r_det * i_B[3] - r_B[3] * i_det ) / mod;
	r_B_inv[1] = -( r_B[1] * r_det + i_B[1] * i_det ) / mod;
	i_B_inv[1] = -( r_det * i_B[1] - r_B[1] * i_det ) / mod;
	r_B_inv[2] = -( r_B[1] * r_det - i_B[1] * i_det ) / mod;
	i_B_inv[2] =  ( r_det * i_B[1] + r_B[1] * i_det ) / mod;
	r_B_inv[3] =  ( r_B[0] * r_det + i_B[0] * i_det ) / mod;
	i_B_inv[3] =  ( r_det * i_B[0] - r_B[0] * i_det ) / mod;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// For generating scale-free graphs: matrix G (p x p) is an adjacency matrix
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void scale_free( int *G, int *p )
{
	int i, j, tmp, dim = *p, p0 = 2;
	double random_value;
	std::vector<int> size_a( dim ); 

	for( i = 0; i < p0 - 1; i++ )
	{
		G[         i * dim + i + 1 ] = 1;
		G[ ( i + 1 ) * dim + i     ] = 1;
	}
		
	for( i = 0 ; i < p0 ; i++ ) size_a[ i ] = 2;
	for( i = p0; i < dim; i++ ) size_a[ i ] = 0;
	
	int total = 2 * p0;
	
	GetRNGstate();
	for( i = p0; i < dim; i++ )
	{
		random_value = (double) total * unif_rand();
	   
		tmp = 0;
		j   = 0;
		
		while( tmp < random_value && j < i ) 
			tmp += size_a[ j++ ];
		
		j--;
		
		G[ i * dim + j ] = 1;
		G[ j * dim + i ] = 1;
		
		total += 2;
		size_a[ j ]++;
		size_a[ i ]++;
	}
	PutRNGstate();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To transfer the raw discreate data 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void transfer_data( int r_data[], int data[], int *n, int *p, int *size_unique_data )
{
	int i, j, l, counter;
	
// - - tranfer each row of raw data to string - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	vector<char> char_row( *p );             
	vector<string>all_patterns( *n );
	string *unique_patterns = new string[ *n ];
	
	for( i = 0; i < *n; i++ )
	{
		for( j = 0; j < *p; j++ )
			char_row[ j ] = r_data[ j * *n + i ] + '0';
		
		all_patterns[ i ] = string( char_row.begin(), char_row.end() );
	}

// - - find the unique string-rows - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	unique_patterns[0] = all_patterns[0];
	int length_unique_patterns = 1;
	for( i = 1; i < *n; i++ )
	{
		counter = 0;
		//for( j = 0; j < length_unique_patterns; j++ )
			//( all_patterns[i] == unique_patterns[j] ) ? j = length_unique_patterns : ++counter;					
		while( ( counter < length_unique_patterns ) and ( all_patterns[ i ] != unique_patterns[ counter ] ) )
			++counter;
		
		if( counter == length_unique_patterns )
			unique_patterns[ length_unique_patterns++ ] = all_patterns[ i ];
	}

// - - tranfer the data - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	int which_one = 0;
	for( l = 0; l < length_unique_patterns; l++ )  
	{
		counter = 0;
		for( i = 0; i < *n; i++ )
			if( all_patterns[ i ] == unique_patterns[ l ] ) 
			{
				counter++;
				which_one = i;
			}
			
		data[ *p * *n + l ] = counter;
		
		for( j = 0; j < *p; j++ )
			data[ j * *n + l ] = r_data[ j * *n + which_one ]; 
	}
	
	*size_unique_data = length_unique_patterns;
	
	delete[] unique_patterns;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

