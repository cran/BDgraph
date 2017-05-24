// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2017 A. (Reza) Mohammadi
//
//     This file is part of BDgraph package.
//
//     BDgraph is free software: you can redistribute it and/or modify it under 
//     the terms of the GNU General Public License as published by the Free 
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
//
//     Maintainer:
//     Reza Mohammadi: a.mohammadi@rug.nl or a.mohammadi@uvt.nl
// ----------------------------------------------------------------------------|
  
#include "matrix.h"

// ----------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves square sub_matrix B (p_sub x p_sub), dictated by vector sub
void sub_matrix( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
	
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
		
		for( j = 0; j < psub; j++ )
			sub_A[ixp + j] = A[subixp + sub[j]]; 
	}
}
   
// ----------------------------------------------------------------------------|
// Takes symmetric matrix A (p x p) and 
// retrieves upper part of sub_matrix B (p_sub x p_sub), dictated by vector sub
void sub_matrix_upper( double A[], double sub_A[], int sub[], int *p_sub, int *p )
{
	int i, j, ixp, subixp, psub = *p_sub, pdim = *p;
			
	for( i = 0; i < psub; i++ )
	{
		ixp    = i * psub;
		subixp = sub[i] * pdim;
			
		for( j = 0; j <= i; j++ )
			sub_A[ixp + j] = A[subixp + sub[j]]; 
	}
}
   
// ----------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves vector sub_A which is 'sub' th row of matrix A, minus 'sub' element
// Likes A[j, -j] in R
void sub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int subj = *sub, pdim = *p, subxp = subj * pdim;

	memcpy( sub_A       , A + subxp           , sizeof( double ) * subj );		
	memcpy( sub_A + subj, A + subxp + subj + 1, sizeof( double ) * ( pdim - subj - 1 ) );	
}
   
// ----------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(2 x p-2) which is sub rows of matrix A, minus two elements
// Likes A[(i,j), -(i,j)] in R ONLY FOR SYMMETRIC MATRICES
void sub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

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

// ----------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves sub_matrix sub_A(p-2 x 2) which is sub cols of matrix A, minus two elements
// Likes A[-(i,j), (i,j)] in R 
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
   
// ----------------------------------------------------------------------------|
// Takes symmatric matrix A (p x p) and 
// retrieves A12(1x(p-1)) and A22((p-1)x(p-1))
// Like A12=A[j, -j], and A22=A[-j, -j] in R
void sub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, ixp1, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;

	memcpy( A12,        A + subxp,            sizeof( double ) * psub );	
	memcpy( A12 + psub, A + subxp + psub + 1, sizeof( double ) * ( pdim - psub - 1 ) );	

	for( i = 0; i < psub; i++ )
	{	
		ixpdim = i * pdim;
		ixp1   = i * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , sizeof( double ) * psub );
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}

	for( i = psub + 1; i < pdim; i++ )
	{
		ixpdim = i * pdim;
		ixp1   = ( i - 1 ) * p1;
		
		memcpy( A22 + ixp1       , A + ixpdim           , sizeof( double ) * psub);
		memcpy( A22 + ixp1 + psub, A + ixpdim + psub + 1, sizeof( double ) * mpsub );
	}
}
    
// ----------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves A11(2x2), A12(2x(p-2)), and A22((p-2)x(p-2))
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
		A12[i + i - 3]       = A[ixp + sub1];
		
		for( j = sub1 + 1; j < pdim; j++ )
			A22[(j - 2) * p2 + i - 2] = A[ixp + j];
	}
}
   
// ----------------------------------------------------------------------------|
// Takes square matrix A (p x p) and 
// retrieves A11_inv(2x2), A21((p-2)x2), and A22((p-2)x(p-2))
// Like A11_inv=inv( A[e, e] ), A21=A[-e, e], and A22=A[-e, -e] in R
void sub_matrices_inv( double A[], double A11_inv[], double A21[], double A22[], int *row, int *col, int *p )
{
	int i, ixp, ixp2, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;
	int sub0xp = sub0 * pdim, sub1xp = sub1 * pdim, sub0_plus = sub0 + 1, sub1_plus = sub1 + 1;
	
	double a11 = A[ sub0 * pdim + sub0 ];
	double a12 = A[ sub0 * pdim + sub1 ];
	double a22 = A[ sub1 * pdim + sub1 ];

	double det_A11 = a11 * a22 - a12 * a12;
	A11_inv[0]     = a22 / det_A11;
	A11_inv[1]     = - a12 / det_A11;
	A11_inv[2]     = A11_inv[1];
	A11_inv[3]     = a11 / det_A11;

	memcpy( A21           , A + sub0xp            , sizeof( double ) * sub0 );		
	memcpy( A21 + sub0    , A + sub0xp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );	
	memcpy( A21 + sub1 - 1, A + sub0xp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	

	memcpy( A21 + p2           , A + sub1xp            , sizeof( double ) * sub0 );		
	memcpy( A21 + p2 + sub0    , A + sub1xp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );	
	memcpy( A21 + p2 + sub1 - 1, A + sub1xp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp  = i * pdim;
		ixp2 = i * p2;

		memcpy( A22 + ixp2           , A + ixp            , sizeof( double ) * sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	
	}
 
	for( i = sub0_plus; i < sub1; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 1 ) * p2;

		memcpy( A22 + ixp2           , A + ixp            , sizeof( double ) * sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );	
	}
	
	for( i = sub1_plus; i < pdim; i++ )
	{
		ixp  = i * pdim;
		ixp2 = ( i - 2 ) * p2;
				
		memcpy( A22 + ixp2           , A + ixp            , sizeof( double ) * sub0 );
		memcpy( A22 + ixp2 + sub0    , A + ixp + sub0_plus, sizeof( double ) * ( sub1 - sub0_plus ) );
		memcpy( A22 + ixp2 + sub1 - 1, A + ixp + sub1_plus, sizeof( double ) * ( pdim - sub1_plus ) );		
	}
}
      
// ----------------------------------------------------------------------------|
// inverse function for symmetric positive-definite matrices (p x p)
// WARNING: Matrix you pass is overwritten with the result
void inverse( double A[], double A_inv[], int *p )
{
	int i, j, info, dim = *p;
	char uplo = 'U';

	// creating an identity matrix
	for( i = 0; i < dim; i++ )
		for( j = 0; j < dim; j++ )
			A_inv[j * dim + i] = (i == j);
	
	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
	F77_NAME(dposv)( &uplo, &dim, &dim, A, &dim, A_inv, &dim, &info );
}
    
// ----------------------------------------------------------------------------|
// inverse function for symmetric (2 x 2)
void inverse_2x2( double B[], double B_inv[] )
{
	double detB = B[0] * B[3] - B[1] * B[1];
	B_inv[0]    = B[3] / detB;
	B_inv[1]    = - B[1] / detB;
	B_inv[2]    = B_inv[1];
	B_inv[3]    = B[0] / detB;
}
    
// ----------------------------------------------------------------------------|
// Cholesky decomposition of symmetric positive-definite matrix
// A = U' %*% U
void cholesky( double A[], double U[], int *p )
{
	char uplo = 'U';
	int j, info, dim = *p, i, pxp = dim * dim;
	
	memcpy( U, A, sizeof( double ) * pxp );	
	
	F77_NAME(dpotrf)( &uplo, &dim, U, &dim, &info );	
	
	for( i = 0; i < dim; i++ )
		for( j = 0; j < i; j++ )
			U[j * dim + i] = 0.0;
}
  
// ----------------------------------------------------------------------------|
// Determinant of symmetric possitive-definite matrix -------------------------|
// ************  WARNING: Matrix you pass is overwritten **********************
//  For any symmetric PD Matrix D, we have:
//                |D| ) = |T| ^ 2
//  where T is the cholesky decomposition of D. Thus, |T| = \prod_{i = 1}^p T_{ii}
//  which makes this quite easy.
void determinant( double A[], double *det_A, int *p )
{
	char uplo = 'U';
	int info, dim = *p, dim1 = dim + 1;
	
	F77_NAME(dpotrf)( &uplo, &dim, &A[0], &dim, &info );	

	double result = 1;
	for( int i = 0; i < dim; i++ ) result *= A[i * dim1];
	
	*det_A = result * result;
}
    
// ----------------------------------------------------------------------------|
// To select an edge for BDMCMC algorithm  
void select_edge( double rates[], int *index_selected_edge, double *sum_rates, int *qp )
{
	int qp_star = *qp;

	// rates = sum_sort_rates
	for ( int i = 1; i < qp_star; i++ )
		rates[i] += rates[ i - 1 ];
	
	*sum_rates   = rates[ qp_star - 1 ];
	double random_value = *sum_rates * runif( 0, 1 );

	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star - 1;
	int position    = upper_bound / 2;  // ( lower_bound + upper_bound ) / 2;

	while( upper_bound - lower_bound > 1 )
	{
		//if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( rates[position] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	*index_selected_edge = ( rates[position] < random_value ) ? ++position : position;
} 
    
// ----------------------------------------------------------------------------|
// To simultaneously select multiple edges for BDMCMC algorithm  
void select_multi_edges( double rates[], int index_selected_edges[], int *size_index, double *sum_rates, int *multi_update, int *qp )
{
	int i, qp_star = *qp, qp_star_1 = qp_star - 1;

	// rates = sum_sort_rates
	for ( i = 1; i < qp_star; i++ )
		rates[i] += rates[ i - 1 ];
	
	double max_bound = rates[qp_star_1];
	
	// ---------- for first edge ----------------------------------------------|
	// To start, find the subscript of the middle position.
	int lower_bound = 0;
	int upper_bound = qp_star_1;
	int position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

	double random_value = max_bound * runif( 0, 1 );

	while( upper_bound - lower_bound > 1 )
	{
		//if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
		( rates[position] > random_value ) ? upper_bound = position : lower_bound = position;     
		
		position = ( lower_bound + upper_bound ) / 2;
	}
	
	if ( rates[position] < random_value ) ++position;
	index_selected_edges[0] = position;
	// ------------------------------------------------------------------------|

	int counter = 1, same;
	for ( int it = 0; it < 200 * *multi_update; it++ )
	{
		if ( counter == *multi_update ) break;
		
		random_value = max_bound * runif( 0, 1 );
	
		// To start, find the subscript of the middle position.
		lower_bound = 0;
		upper_bound = qp_star_1;
		position    = upper_bound / 2; // ( lower_bound + upper_bound ) / 2;

		while( upper_bound - lower_bound > 1 )
		{
			//if ( rates[position] > random_value ) { upper_bound = position; } else { lower_bound = position; }     
			( rates[position] > random_value ) ? upper_bound = position : lower_bound = position;     
			
			position = ( lower_bound + upper_bound ) / 2;
		}
		
		if ( rates[position] < random_value ) ++position;
		
		same = 0;
		for ( i = 0; i < counter; i++ )
			if( index_selected_edges[i] == position )
				++same;

		if ( same == 0 ) index_selected_edges[counter++] = position;
	}
	
	*size_index = counter;
	*sum_rates  = max_bound;
} 
         
// ----------------------------------------------------------------------------|
// Parallel Computation for birth-death rates for BD-MCMC algorithm
// ----------------------------------------------------------------------------|
void rates_bdmcmc_parallel( double rates[], int G[], int index_row[], int index_col[], int *sub_qp, double Ds[], double Dsijj[],
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
			i = index_row[ counter ];
			j = index_col[ counter ];

			jj   = j * dim1;
			Dsjj = Ds[jj];
			
			sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], &j, &dim );

			// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
			// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
			sigmajj_inv = - 1.0 / sigma[jj];
			F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 );

			ij = j * dim + i;
			
			// For (i,j) = 0 ----------------------------------------------|	
			sub_row_mins( &K[0], &Kj12[0], &j, &dim );   // Kj12 = K[j, -j]  
			Kj12[ i ] = 0.0;                         // Kj12[1,i] = 0

			// Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
			F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one );
			
			// K022 = Kj12xK22_inv %*% t(Kj12)
			K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			

			// For (i,j) = 1 ----------------------------------------------|
			sub_cols_mins( &K[0], &K21[0], &i, &j, &dim );  // K21 = K[-e, e]  
			
			sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], &i, &j, &dim );

			// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
			F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 );

			// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
			F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 );

			// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
			F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two );  
			
			// K121 = K12xK22_inv %*% K21													
			F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two );		
			// Finished (i,j) = 1------------------------------------------|

			a11      = K[i * dim1] - K121[0];	
			sum_diag = Dsjj * ( K022 - K121[3] ) - Ds[ij] * ( K121[1] + K121[2] );

			// nu_star = b + sum( Gf[,i] * Gf[,j] )
			nu_star = b1;
			//for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k];   
			for( k = 0; k < dim; k++ ) nu_star += G[i * dim + k] * G[j * dim + k];   
			//nu_star = F77_NAME(ddot)( &dim, &G[0] + ixdim, &one, &G[0] + jxdim, &one );
			nu_star = 0.5 * nu_star;

			log_rate = ( G[ij] )   
				? 0.5 * log( 2.0 * Dsjj / a11 ) + lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsijj[ij] * a11 + sum_diag )
				: 0.5 * log( 0.5 * a11 / Dsjj ) - lgammafn( nu_star + 0.5 ) + lgammafn( nu_star ) + 0.5 * ( Dsijj[ij] * a11 + sum_diag );

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
     
// ----------------------------------------------------------------------------|
// computing birth/death rate or alpha for element (i,j)
// it is for double Metropolis-Hasting algorihtms
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
	F77_NAME(dsyr)( &sideL, p1, &sigmajj_inv, sigmaj12, &one, sigmaj22, p1 );

	// For (i,j) = 0 ----------------------------------------------|	
	sub_row_mins( K, Kj12, selected_edge_j, dim );  // K12 = K[j, -j]  
	Kj12[ *selected_edge_i ] = 0.0;                       // K12[1,i] = 0

	// Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
	F77_NAME(dsymv)( &sideL, p1, &alpha, &sigmaj22[0], p1, Kj12, &one, &beta, Kj12xK22_inv, &one );
	
	// K022 = Kj12xK22_inv %*% t(Kj12)
	double K022 = F77_NAME(ddot)( p1, Kj12xK22_inv, &one, Kj12, &one );			

	// For (i,j) = 1 ----------------------------------------------|
	sub_cols_mins( K, K12, selected_edge_i, selected_edge_j, dim );   // K21 = K[-e, e] 
	
	sub_matrices_inv( sigma, sigma11_inv, sigma12, sigma22, selected_edge_i, selected_edge_j, dim );

	// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
	F77_NAME(dgemm)( &transN, &transN, p2, &two, &two, &alpha, sigma12, p2, sigma11_inv, &two, &beta, sigma21xsigma11_inv, p2 );

	// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
	F77_NAME(dgemm)( &transN, &transT, p2, p2, &two, &alpha1, sigma21xsigma11_inv, p2, sigma12, p2, &beta1, sigma22, p2 );

	// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
	F77_NAME(dgemm)( &transT, &transN, &two, p2, p2, &alpha, K12, p2, sigma22, p2, &beta, K12xK22_inv, &two );  
	
	// K121 = K12xK22_inv %*% K21													
	F77_NAME(dgemm)( &transN, &transN, &two, &two, p2, &alpha, K12xK22_inv, &two, K12, p2, &beta, K121, &two );		
	// Finished (i,j) = 1------------------------------------------|

	double a11      = K[*selected_edge_i * *dim + *selected_edge_i] - K121[0];	
	double sum_diag = *Dsjj * ( K022 - K121[3] ) - *Dsij * ( K121[1] + K121[2] );

	// Dsijj = Dsii - Dsij * Dsij / Dsjj;
	*log_Hij = ( log( static_cast<double>(*Dsjj) ) - log( static_cast<double>(a11) ) + *Dsijj * a11 - sum_diag ) / 2;
}    
     
// ----------------------------------------------------------------------------|
// Parallel Computation for birth-death rates for double BD-MCMC algorithm
// ----------------------------------------------------------------------------|
void rates_bdmcmc_dmh_parallel( double rates[], int G[], int index_row[], int index_col[], int *sub_qp, double Ds[], double D[],
				            double sigma[], double K[], double sigma_dmh[], 
				            double K_dmh[], int *b, int *p )
{
	int dim = *p, p1 = dim - 1, p2 = dim - 2, p2x2 = ( dim - 2 ) * 2;

	#pragma omp parallel
	{
		int index_rate_j, i, j, ij, jj;
		double Dsjj, Dsij, Dsijj, Dij, Dijj, Djj, log_rate;

		double *K121         = new double[ 4 ];  
		double *Kj12         = new double[ p1 ];  
		double *sigmaj12     = new double[ p1 ];  
		double *sigmaj22     = new double[ p1 * p1 ];  
		double *Kj12xK22_inv = new double[ p1 ];  
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
			Dsjj = Ds[jj];
			Djj  = D[jj];

			for( i = 0; i < j; i++ )
			{
				ij    = j * dim + i;
				Dsij  = Ds[ij];
				Dsijj = - Dsij * Dsij / Dsjj;
				Dij   = D[ij];
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
				
				log_rate = ( G[ij] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );				
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
     	
// -------------- NEW for Lang codes ------------------------------------------|
// For Hermitian matrix
void Hsub_row_mins( double A[], double sub_A[], int *sub, int *p )
{
	int i, l = 0, subj = *sub, pdim = *p, subxp = subj * pdim;

	for( i = 0; i < subj; i++ )
		sub_A[l++] = -A[subxp + i];
	
	for( i = subj + 1; i < pdim; i++ )
		sub_A[l++] = -A[subxp + i];
}
      
// ----------------------------------------------------------------------------|
// For Hermitian matrix
void Hsub_rows_mins( double A[], double sub_A[], int *row, int *col, int *p )
{	
	int i, l = 0, pdim = *p, sub0 = *row, sub1 = *col, sub0p = sub0 * pdim, sub1p = sub1 * pdim;

	for( i = 0; i < sub0; i++ )
	{
		sub_A[l++] = -A[sub0p + i]; 
		sub_A[l++] = -A[sub1p + i]; 
	}
	
	for( i = sub0 + 1; i < sub1; i++ )
	{
		sub_A[l++] = -A[sub0p + i]; 
		sub_A[l++] = -A[sub1p + i]; 
	}

	for( i = sub1 + 1; i < pdim; i++ )
	{
		sub_A[l++] = -A[sub0p + i]; 
		sub_A[l++] = -A[sub1p + i]; 
	}
}
       
// ----------------------------------------------------------------------------|
// sub_matrices1 for Hermitian matrix
void Hsub_matrices1( double A[], double A12[], double A22[], int *sub, int *p )
{
	int i, ixpdim, pdim = *p, p1 = pdim - 1, psub = *sub, subxp = psub * pdim, mpsub = pdim - psub - 1;

	for( i = 0; i < psub; i++ )
		A12[i] = -A[subxp + i];
	for( i = psub; i < pdim - 1; i++ )
		A12[i] = -A[subxp + i + 1];

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
        
// ----------------------------------------------------------------------------|
// sub_matrices for Hermitian matrix
void Hsub_matrices( double A[], double A11[], double A12[], double A22[], int *row, int *col, int *p )
{
	int i, i1, i2, ixp, pdim = *p, p2 = pdim - 2, sub0 = *row, sub1 = *col;

	A11[0] = A[sub0 * pdim + sub0];
	A11[1] = A[sub0 * pdim + sub1];
	A11[2] = -A11[1];                   // for symmetric matrices
	A11[3] = A[sub1 * pdim + sub1];
 
	for( i = 0; i < sub0; i++ )
	{	
		ixp = i * pdim;
		
		A12[i + i]     = A[ixp + sub0];
		A12[i + i + 1] = A[ixp + sub1];

		memcpy( A22 + i * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
 
	for( i = sub0 + 1; i < sub1; i++ )
	{
		ixp = i * pdim;
		i1 = i - 1;

		A12[i + i - 2] = A[ixp + sub0];
		A12[i + i - 1] = A[ixp + sub1];

		memcpy( A22 + i1 * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i1 * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i1 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );	
	}
	
	for( i = sub1 + 1; i < pdim; i++ )
	{
		ixp = i * pdim;
		i2  = i - 2;
				
		A12[i + i - 4] = A[ixp + sub0];
		A12[i + i - 3] = A[ixp + sub1];

		memcpy( A22 + i2 * p2,            A + ixp,            sizeof( double ) * sub0 );
		memcpy( A22 + i2 * p2 + sub0,     A + ixp + sub0 + 1, sizeof( double ) * ( sub1 - sub0 - 1 ) );
		memcpy( A22 + i2 * p2 + sub1 - 1, A + ixp + sub1 + 1, sizeof( double ) * ( pdim - sub1 - 1 ) );		
	}
}
   
// ----------------------------------------------------------------------------|
// inverse function for Hermitian (2 x 2)
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

// ----------------------------------------------------------------------------|
// For generating scale-free graphs: matrix G (p x p) is an adjacency matrix
void scale_free( int *G, int *p )
{
	GetRNGstate();
	int i, j, tmp, total, dim = *p, p0 = 2;
	double random_value;
	std::vector<int> size_a( dim ); 

	for( i = 0; i < p0 - 1; i++ )
	{
		G[i * dim + i + 1]   = 1;
		G[(i + 1) * dim + i] = 1;
	}
		
	for( i = 0; i < p0; i++ ) size_a[i] = 2;
	
	for( i = p0; i < dim; i++ ) size_a[i] = 0;
	
	total = 2 * p0;
	
	for( i = p0; i < dim; i++ )
	{
	   random_value = (double) total * runif( 0, 1 );
	   
		tmp = 0;
		j   = 0;
		
		while( tmp < random_value && j < i ) 
			tmp += size_a[j++];
		
		j--;
		G[i * dim + j] = 1;
		G[j * dim + i] = 1;
		total += 2;
		size_a[j]++;
		size_a[i]++;
	}
	PutRNGstate();
}
