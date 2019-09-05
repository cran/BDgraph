// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2019  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  
#include "copula.h"
   
extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing mean for copula function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_mean( double Z[], double K[], double *mu_ij, double *sigma, int *i, int *j, int *n, int *p )
{
    int k, dim = *p, number = *n, row = *i, col = *j, jxp = col * dim;
    double mu = 0.0;
    
    for( k = 0;       k < col; k++ ) mu += Z[ k * number + row ] * K[ jxp + k ];
    for( k = col + 1; k < dim; k++ ) mu += Z[ k * number + row ] * K[ jxp + k ];
    
    *mu_ij = - mu * *sigma;
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing bounds for copula function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_bounds( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
    int kj, col = *j, number = *n, ij = col * number + *i;
    double low_b = -1e308, upper_b = +1e308;
    
    for( int k = 0; k < number; k++ )
    {
        kj = col * number + k;
        
        if( R[ kj ] < R[ ij ] ) 
            low_b = max( Z[ kj ], low_b );	
        else if( R[ kj ] > R[ ij ] ) 
            upper_b = min( Z[ kj ], upper_b );
    }
    
    *lb = low_b;
    *ub = upper_b;	
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// copula for BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void copula( double Z[], double K[], int R[], int not_continuous[], int *n, int *p )
{
    int number = *n, dim = *p, nxp = number * dim, dimp1 = dim + 1;
    
    #pragma omp parallel
    {	
        double sigma, sd_j, mu_ij, lb, ub, runif_value, pnorm_lb, pnorm_ub;
        int i, j;
        
        #pragma omp for
        for( int counter = 0; counter < nxp; counter++ )
        {   
            j = counter / number;
            i = counter % number;
            
            if( not_continuous[ j ] )
            {
                sigma = 1.0 / K[ j * dimp1 ]; // 1.0 / K[ j * dim + j ];
                sd_j  = sqrt( sigma );
                
                get_mean( Z, K, &mu_ij, &sigma, &i, &j, &number, &dim );
                
                get_bounds_NA( Z, R, &lb, &ub, &i, &j, &number );
                
                pnorm_lb     = Rf_pnorm5( lb, mu_ij, sd_j, TRUE, FALSE );
                pnorm_ub     = Rf_pnorm5( ub, mu_ij, sd_j, TRUE, FALSE );
                //runif_value = runif( pnorm_lb, pnorm_ub );
                runif_value  = pnorm_lb + unif_rand() * ( pnorm_ub - pnorm_lb );
                Z[ counter ] = Rf_qnorm5( runif_value, mu_ij, sd_j, TRUE, FALSE );
            }
        }
    }
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing bounds for copula function for data with missing values 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_bounds_NA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
    int kj, number = *n, col = *j, ij = col * number + *i;
    double low_b = -1e308, upper_b = +1e308;
    
    for( int k = 0; k < number; k++ )
    {
        kj = col * number + k;
       
        if( R[ kj ] != 0 )
        {
            if( R[ kj ] < R[ ij ] ) 
                low_b = max( Z[ kj ], low_b );	
            else if( R[ kj ] > R[ ij ] ) 
                upper_b = min( Z[ kj ], upper_b );
        }
    }
    
    *lb = low_b;
    *ub = upper_b;		
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// copula for data with missing values 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void copula_NA( double Z[], double K[], int R[], int not_continuous[], int *n, int *p )
{
    int number = *n, dim = *p, nxp = number * dim, dimp1 = dim + 1;
    
    #pragma omp parallel
    {	
        double sigma, sd_j, mu_ij, lb, ub, runif_value, pnorm_lb, pnorm_ub;
        int i, j;
        
        #pragma omp for
        for( int counter = 0; counter < nxp; counter++ )
        {   
            j = counter / number;
            i = counter % number;
            
            if( not_continuous[ j ] )
            {
                sigma = 1.0 / K[ j * dimp1 ]; // 1.0 / K[ j * dim + j ];
                sd_j  = sqrt( sigma );
                
                get_mean( Z, K, &mu_ij, &sigma, &i, &j, &number, &dim );
                
                if( R[ counter ] != 0 )
                {
                    get_bounds_NA( Z, R, &lb, &ub, &i, &j, &number );
                    
                    pnorm_lb     = Rf_pnorm5( lb, mu_ij, sd_j, TRUE, FALSE );
                    pnorm_ub     = Rf_pnorm5( ub, mu_ij, sd_j, TRUE, FALSE );
                    //runif_value = runif( pnorm_lb, pnorm_ub );
                    runif_value  = pnorm_lb + unif_rand() * ( pnorm_ub - pnorm_lb );
                    Z[ counter ] = Rf_qnorm5( runif_value, mu_ij, sd_j, TRUE, FALSE );
                }else
                    Z[ counter ] = mu_ij + norm_rand() * sd_j;  // rnorm( mu_ij, sd_j );
            }
        }
    }
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Calculating Ds = D + S for the BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_Ds( double K[], double Z[], int R[], int not_continuous[], double D[], double Ds[], double S[], int *gcgm, int *n, int *p )
{
	int dim = *p;

	( *gcgm == 0 ) ? copula( Z, K, R, not_continuous, n, &dim ) : copula_NA( Z, K, R, not_continuous, n, &dim );
	
	// S <- t(Z) %*% Z; NOTE, I'm using Ds instead of S, for saving memory
	double alpha = 1.0, beta  = 0.0;
	char transA = 'T', transB = 'N';
	F77_NAME(dgemm)( &transA, &transB, &dim, &dim, n, &alpha, Z, n, Z, n, &beta, &S[0], &dim FCONE FCONE );		
	
    #pragma omp parallel for
	for( int i = 0; i < dim * dim; i++ ) 
	    Ds[ i ] = D[ i ] + S[ i ];		
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Calculating Ts = chol( solve( Ds ) ) for the BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_Ts( double Ds[], double Ts[], double inv_Ds[], double copy_Ds[], int *p )
{
	int dim = *p; 

	memcpy( &copy_Ds[0], Ds, sizeof( double ) * dim * dim );
	
	inverse( &copy_Ds[0], &inv_Ds[0], &dim );	

	cholesky( &inv_Ds[0], Ts, &dim );	
}

} // End of exturn "C"
