// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2016 Mohammadi A. and Wit C. E.
//
//     This file is part of BDgraph package.
//
//     BDgraph is free software: you can redistribute it and/or modify it under 
//     the terms of the GNU General Public License as published by the Free 
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
//
//     Maintainer:
//     Abdolreza Mohammadi: a.mohammadi@rug.nl or a.mohammadi@uvt.nl
// ----------------------------------------------------------------------------|
  
#include "copula.h"
   
// Calculating mean for copula function
void get_mean( double Z[], double K[], double *mu_ij, double *sigma, int *i, int *j, int *n, int *p )
{
	int k, dim = *p, number = *n, row = *i, col = *j;
	double mu = 0.0;
	
	for( k = 0; k < col; k++ ) 
		mu += Z[k * number + row] * K[col * dim + k];

	for( k = col + 1; k < dim; k++ ) 
		mu += Z[k * number + row] * K[col * dim + k];

	*mu_ij = - mu * *sigma;
}

// Calculating bounds for copula function
void get_bounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
	int kj, ij, row = *i, col = *j;
	double low_b = -1e308, upper_b = +1e308;

	for( int k = 0, number = *n; k < number; k++ )
	{
		kj = col * number + k;
		ij = col * number + row;
		
     // if( R[k, j] < R[i, j] ) lb = max( Z[ k, j], lb )
	 // if( R[k, j] > R[i, j] ) ub = min( Z[ k, j], ub )										
		if( R[kj] < R[ij] ) 
			low_b = max( Z[kj], low_b );	
		else if( R[kj] > R[ij] ) 
			upper_b = min( Z[kj], upper_b );
	}

	*lb = low_b;
	*ub = upper_b;	
}
 
// copula part
void copula( double Z[], double K[], int R[], int *n, int *p )
{
	int number = *n, dim = *p;
	double sigma, sd_j, mu_ij;
	
	double lb, ub;
	double runif_value, pnorm_lb, pnorm_ub;
	
	for( int j = 0; j < dim; j++ )
	{   
		sigma = 1 / K[j * dim + j];
		sd_j   = sqrt( sigma );
		
		// interval and sampling
		for( int i = 0; i < number; i++ )
		{
			get_mean( Z, K, &mu_ij, &sigma, &i, &j, &number, &dim );
			
			get_bounds( Z, R, &lb, &ub, &i, &j, &number );
			
			// runif_value = runif( 1, pnorm( lower_bound, mu_ij, sd_j ), pnorm( upper_bound, mu_ij, sd_j ) )
			// Z[i,j]     = qnorm( runif_value, mu_ij, sd_j )									
			//GetRNGstate();
			pnorm_lb          = pnorm( lb, mu_ij, sd_j, TRUE, FALSE );
			pnorm_ub          = pnorm( ub, mu_ij, sd_j, TRUE, FALSE );
			runif_value       = runif( pnorm_lb, pnorm_ub );
			Z[j * number + i] = qnorm( runif_value, mu_ij, sd_j, TRUE, FALSE );
			//PutRNGstate();				
		}
	}
}

// Calculating bounds for copula function with missing data 
void get_bounds_NA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
	int kj, ij, row = *i, col = *j;
	double low_b = -1e308, upper_b = +1e308;

	for( int k = 0, number = *n; k < number; k++ )
	{
		kj = col * number + k;
		ij = col * number + row;
		
		if( R[kj] != 0 )
		{
			// if( R[k, j] < R[i, j] ) lb = max( Z[ k, j], lb )
			// if( R[k, j] > R[i, j] ) ub = min( Z[ k, j], ub )										
			if( R[kj] < R[ij] ) 
				low_b = max( Z[kj], low_b );	
			else if( R[kj] > R[ij] ) 
				upper_b = min( Z[kj], upper_b );
		}
	}
	
	*lb = low_b;
	*ub = upper_b;		
}
 
// copula part for missing data
void copula_NA( double Z[], double K[], int R[], int *n, int *p )
{
	int ij, number = *n, dim = *p;
	double sigma, sd_j, mu_ij, lb, ub, runif_value, pnorm_lb, pnorm_ub;
	
	for( int j = 0; j < dim; j++ )
	{   
		sigma = 1 / K[j * dim + j];
		sd_j   = sqrt( sigma );
		
		// interval and sampling
		for( int i = 0; i < number; i++ )
		{
			get_mean( Z, K, &mu_ij, &sigma, &i, &j, &number, &dim );
			
			ij = j * number + i;
			//GetRNGstate();
			if( R[ij] != 0 )
			{
				get_bounds_NA( Z, R, &lb, &ub, &i, &j, &number );
				
				pnorm_lb    = pnorm( lb, mu_ij, sd_j, TRUE, FALSE );
				pnorm_ub    = pnorm( ub, mu_ij, sd_j, TRUE, FALSE );
				runif_value = runif( pnorm_lb, pnorm_ub );
				Z[ij]       = qnorm( runif_value, mu_ij, sd_j, TRUE, FALSE );
			} 
			else 
				Z[ij] = rnorm( mu_ij, sd_j );
			//PutRNGstate();				
		}
	}
}
    
// Calculating Ds = D + S for the BDMCMC sampling algorithm
void get_Ds( double K[], double Z[], int R[], double D[], double Ds[], double S[], int *gcgm, int *n, int *p )
{
	int gcgm_check = *gcgm, dim = *p, pxp = dim * dim;

	if( gcgm_check == 0 )
		copula( Z, K, R, n, &dim );
	else
		copula_NA( Z, K, R, n, &dim );
	
	// S <- t(Z) %*% Z
	// Here, I'm using Ds instead of S, for saving memory
	double alpha = 1.0, beta  = 0.0;
	char transA = 'T', transB = 'N';
	// LAPACK function to compute  C := alpha * A * B + beta * C
	//        DGEMM ( TRANSA,  TRANSB, M, N, K,  ALPHA, A,LDA,B, LDB,BETA, C, LDC )																				
	F77_NAME(dgemm)( &transA, &transB, &dim, &dim, n, &alpha, Z, n, Z, n, &beta, &S[0], &dim );		
	// Ds = D + S
	for( int i = 0; i < pxp ; i++ ) Ds[i] = D[i] + S[i];		
}

// Calculating Ts = chol( solve( Ds ) ) for the BDMCMC sampling algorithm
void get_Ts( double Ds[], double Ts[], double inv_Ds[], double copy_Ds[], int *p )
{
	int dim = *p, pxp = dim * dim;

	memcpy( &copy_Ds[0], Ds, sizeof( double ) * pxp );
	
	inverse( &copy_Ds[0], &inv_Ds[0], &dim );	

	cholesky( &inv_Ds[0], Ts, &dim );	
}
