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
  
#include "copula.h"
   
extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing mean for copula function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_mean( double Z[], double K[], double *mu_ij, double *sigma, int *i, int *j, int *n, int *p )
{
    int k, dim = *p, number = *n, row = *i, col = *j, jxp = col * dim;
    double mu = 0.0;
    
    for( k = 0;       k < col; k++ ) mu += Z[ k * number + row ] * K[ jxp + k ];
    for( k = col + 1; k < dim; k++ ) mu += Z[ k * number + row ] * K[ jxp + k ];
    
    *mu_ij = - mu * *sigma;
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing bounds for copula function
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
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
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// copula for BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
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
                
                get_bounds( Z, R, &lb, &ub, &i, &j, &number );
                
                pnorm_lb     = Rf_pnorm5( lb, mu_ij, sd_j, TRUE, FALSE );
                pnorm_ub     = Rf_pnorm5( ub, mu_ij, sd_j, TRUE, FALSE );
                //runif_value = runif( pnorm_lb, pnorm_ub );
                runif_value  = pnorm_lb + unif_rand() * ( pnorm_ub - pnorm_lb );
                Z[ counter ] = Rf_qnorm5( runif_value, mu_ij, sd_j, TRUE, FALSE );
            }
        }
    }
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
// copula - Discrete Weibull for BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
void copula_dw( double Z[], double K[], int Y[], double lower_bounds[], double upper_bounds[], 
                int *n, int *p )
{
    int number = *n, dim = *p, nxp = number * dim, dimp1 = dim + 1;
    
    #pragma omp parallel
    {	
        double sigma, sd_j, mu_ij, runif_value, pnorm_lb, pnorm_ub;
        int i, j;
        
        #pragma omp for
        for( int counter = 0; counter < nxp; counter++ )
        {   
            j = counter / number;
            i = counter % number;
            
            sigma = 1.0 / K[ j * dimp1 ]; // 1.0 / K[ j * dim + j ];
            sd_j  = sqrt( sigma );
            
            get_mean( Z, K, &mu_ij, &sigma, &i, &j, &number, &dim );
            
            pnorm_lb     = Rf_pnorm5( lower_bounds[ counter ], mu_ij, sd_j, TRUE, FALSE );
            pnorm_ub     = Rf_pnorm5( upper_bounds[ counter ], mu_ij, sd_j, TRUE, FALSE );
            //runif_value = runif( pnorm_lb, pnorm_ub );
            runif_value  = pnorm_lb + unif_rand() * ( pnorm_ub - pnorm_lb );
            Z[ counter ] = Rf_qnorm5( runif_value, mu_ij, sd_j, TRUE, FALSE );
        }
    }
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
// copula - Discrete Weibull for data with missing values 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
void copula_dw_NA( double Z[], double K[], int Y[], double lower_bounds[], double upper_bounds[], 
                   int *n, int *p )
{
    int number = *n, dim = *p, nxp = number * dim, dimp1 = dim + 1;
    
    #pragma omp parallel
    {	
        double sigma, sd_j, mu_ij, runif_value, pnorm_lb, pnorm_ub;
        int i, j;
        
        #pragma omp for
        for( int counter = 0; counter < nxp; counter++ )
        {   
            j = counter / number;
            i = counter % number;
            
            sigma = 1.0 / K[ j * dimp1 ]; // 1.0 / K[ j * dim + j ];
            sd_j  = sqrt( sigma );
            
            get_mean( Z, K, &mu_ij, &sigma, &i, &j, &number, &dim );
            
            if( Y[ counter ] != -1000 ) // here NA values have been replaced by -1000
            {
                pnorm_lb     = Rf_pnorm5( lower_bounds[ counter ], mu_ij, sd_j, TRUE, FALSE );
                pnorm_ub     = Rf_pnorm5( upper_bounds[ counter ], mu_ij, sd_j, TRUE, FALSE );
                //runif_value = runif( pnorm_lb, pnorm_ub );
                runif_value  = pnorm_lb + unif_rand() * ( pnorm_ub - pnorm_lb );
                Z[ counter ] = Rf_qnorm5( runif_value, mu_ij, sd_j, TRUE, FALSE );
            }else
                Z[ counter ] = mu_ij + norm_rand() * sd_j;  // rnorm( mu_ij, sd_j );
        }
    }
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Calculating Ds = D + S for the copula-Discrete Weibull in BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_Ds_dw( double K[], double Z[], int Y[], double lower_bounds[], 
                double upper_bounds[], double D[], double Ds[], double S[], 
                int *gcgm, int *n, int *p )
{
    int dim = *p;
    
    ( *gcgm == 0 ) ? copula_dw( Z, K, Y, lower_bounds, upper_bounds, n, &dim ) : copula_dw_NA( Z, K, Y, lower_bounds, upper_bounds, n, &dim );
    
    // S <- t(Z) %*% Z; NOTE, I'm using Ds instead of S, for saving memory
    double alpha = 1.0, beta  = 0.0;
    char transA = 'T', transB = 'N';
    F77_NAME(dgemm)( &transA, &transB, &dim, &dim, n, &alpha, Z, n, Z, n, &beta, &S[0], &dim FCONE FCONE );		
    
    #pragma omp parallel for
    for( int i = 0; i < dim * dim; i++ ) 
        Ds[ i ] = D[ i ] + S[ i ];		
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing bounds for copula function for data with missing values 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_bounds_NA( double Z[], int R[], double *lb, double *ub, int *i, int *j, int *n )
{
    int kj, number = *n, col = *j, ij = col * number + *i;
    double low_b = -1e308, upper_b = +1e308;
    
    for( int k = 0; k < number; k++ )
    {
        kj = col * number + k;
       
        if( R[ kj ] != -1000 )  // here NA values have been replaced by -1000 
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
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// copula for data with missing values 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
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
                
                if( R[ counter ] != -1000 ) // here NA values have been replaced by -1000
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
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Calculating Ds = D + S for the BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_Ds( double K[], double Z[], int R[], int not_continuous[], double D[], 
             double Ds[], double S[], int *gcgm, int *n, int *p )
{
	int dim = *p;

	( *gcgm == 0 ) ? copula( Z, K, R, not_continuous, n, &dim ) : copula_NA( Z, K, R, not_continuous, n, &dim );
	
	// S <- t(Z) %*% Z; NOTE, I use Ds instead of S, to save memory
	double alpha = 1.0, beta  = 0.0;
	char transA = 'T', transB = 'N';
	F77_NAME(dgemm)( &transA, &transB, &dim, &dim, n, &alpha, Z, n, Z, n, &beta, &S[0], &dim FCONE FCONE );		
	
    #pragma omp parallel for
	for( int i = 0; i < dim * dim; i++ ) 
	    Ds[ i ] = D[ i ] + S[ i ];		
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Calculating Ts = chol( solve( Ds ) ) for the BDMCMC sampling algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_Ts( double Ds[], double Ts[], double inv_Ds[], double copy_Ds[], int *p )
{
	int dim = *p; 

	memcpy( &copy_Ds[0], Ds, sizeof( double ) * dim * dim );
	
	inverse( &copy_Ds[0], &inv_Ds[0], &dim );	

	cholesky( &inv_Ds[0], Ts, &dim );	
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To update tu for the tgm models
//	d_mu_i    = data[ i, , drop = FALSE ] - mu       # 1 x p
//	delta_y_i = c( d_mu_i %*% K %*% t( d_mu_i ) )    # 1 x 1
	
//	shape_tu_i = ( nu + p ) / 2
//	rate_tu_i  = ( nu + delta_y_i ) / 2
	
//	tu[ i ] = rgamma( 1, shape = shape_tu_i, rate = rate_tu_i )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_tu( double data[], double K[], double tu[], double mu[], double *nu, int *n, int *p )
{
    int i, j, k, l, dim = *p, size_data = *n; //, one = 1;
    double nu_c = *nu, delta_y_i, shape_tu_i, rate_tu_i;
    
	vector<double> d_mu_i( dim ); 
	//vector<double> d_mu_i_x_K( dim ); 

    //GetRNGstate();
    for( i = 0; i < size_data; i++ )
    {
        // for( j in 1:p ) d_mu_i[ j ] = data[ i, j ] - mu[ j ];
        for( j = 0; j < dim; j++ )
             d_mu_i[ j ] = data[ j * size_data + i ] - mu[ j ];
      
        // delta_y_i = 0.0;
        // for( k in 1:p )
        //    for( l in 1:p )
        //        delta_y_i = delta_y_i + d_mu_i[ l ] * K[ l, k ] * d_mu_i[ k ];
      
        delta_y_i = 0.0;
        for( k = 0; k < dim; k++ )
            for( l = 0; l < dim; l++ )
                delta_y_i += d_mu_i[ l ] * K[ k * dim + l ] * d_mu_i[ k ];
        
        shape_tu_i = ( nu_c + static_cast<double>( dim ) ) / 2.0;
        rate_tu_i  = ( nu_c + delta_y_i ) / 2.0;
                
        // tu[ i ]   = rgamma( 1, shape = shape_tu_i, scale = 1.0 / rate_tu_i )
        tu[ i ] = Rf_rgamma( shape_tu_i, 1.0 / rate_tu_i );
    }
    //PutRNGstate();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To update Ds for the tgm models
// Ds = D + S for the BDMCMC sampling algorithm
//	for( i in 1:p )
//    for( j in 1:p )
//        for( k in 1:n )
//            S[ i, j ] = S[ i, j ] + tu[ k ] * ( data[ k, i ] - mu[ i ] ) * ( data[ k, j ] - mu[ j ] )
//    	Ds = D + S
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void get_Ds_tgm( double data[], double D[], double mu[], double tu[], double Ds[], double S[], int *n, int *p )
{
	int i, j, k, ij, dim = *p, size_data = *n;

    for( i = 0; i < dim; i++ )
        for( j = 0; j < dim; j++ )
        {
            ij = j * dim + i;
            for( k = 0; k < size_data; k++ )
                S[ ij ] += tu[ k ] * ( data[ i * size_data + k ] - mu[ i ] ) * ( data[ j * size_data + k ] - mu[ j ] );
        }

    #pragma omp parallel for
	for( int i = 0; i < dim * dim; i++ ) 
	    Ds[ i ] = D[ i ] + S[ i ];		
}
   
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// To update mu for the tgm models
// mu = tu %*% data / sum( tu )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void update_mu( double data[], double mu[], double tu[], int *n, int *p )
{
   	int i, j, dim = *p, one = 1;
    
   	// mu = tu %*% data / sum( tu )
    double alpha = 1.0, beta  = 0.0;
    char transA = 'N', transB = 'N';
    //F77_NAME(dgemm)( &transA, &transB, &one, &dim, &size_data, &alpha, &tu_c[0], &one, &data_c[0], &size_data, &beta, &mu_c[0], &one FCONE FCONE );		
    F77_NAME(dgemm)( &transA, &transB, &one, &dim, n, &alpha, &tu[0], &one, &data[0], n, &beta, &mu[0], &one FCONE FCONE );		

    double sum_tu = 0.0;
    for( i = 0; i < *n; i++ )
       sum_tu += tu[ i ];
	
	// mu = mu / sum( tu )
    for( j = 0; j < dim; j++ )
        mu[ j ] = mu[ j ] / sum_tu;
}
  
} // End of exturn "C"
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
