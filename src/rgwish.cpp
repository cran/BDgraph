// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2020  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  
#include "rgwish.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// sampling from Wishart distribution, in which Ts = chol( solve( Ds ) )
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rwish_c( double Ts[], double K[], int *b, int *p )
{
	int dim = *p, pxp = dim * dim, bK = *b;
	double alpha = 1.0, beta   = 0.0;
	char transT  = 'T', transN = 'N', side = 'R', upper = 'U';																	

	vector<double> psi( pxp, 0.0 ); 

	// - - - Sample values in Psi matrix - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	GetRNGstate();
	//#pragma omp parallel for
	for( int i = 0; i < dim; i++ )
		psi[ i * dim + i ] = sqrt( Rf_rgamma( ( bK + dim - i - 1 ) / 2.0, 2.0 ) );
		//psi[i * dim + i] = sqrt( rchisq( bK + dim - i - 1 ) );

	//#pragma omp parallel for
	for( int j = 1; j < dim; j++ )
		for( int i = 0; i < j; i++ )
			psi[ j * dim + i ] = norm_rand();
	PutRNGstate();
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

    	// C = psi %*% Ts   I used   psi = psi %*% Ts
	// dtrmm (SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
	F77_NAME(dtrmm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &psi[0], &dim FCONE FCONE FCONE FCONE );

	// K = t(C) %*% C 
	// LAPACK function to compute  C := alpha * A * B + beta * C																				
	F77_NAME(dgemm)( &transT, &transN, &dim, &dim, &dim, &alpha, &psi[0], &dim, &psi[0], &dim, &beta, K, &dim FCONE FCONE );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// G is adjacency matrix which has zero in its diagonal // threshold = 1e-8
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rgwish_c( int G[], double Ts[], double K[], int *b, int *p, double *threshold )
{
	int info, i, j, l, size_node, one = 1, dim = *p, pxp = dim * dim;	
    
    double threshold_c = *threshold;
	double alpha = 1.0, beta  = 0.0;
	
	char transN  = 'N', uplo  = 'U'; 
	
	rwish_c( Ts, K, b, &dim );
	
	vector<double> sigma_start( pxp ); 
	inverse( K, &sigma_start[0], &dim );
	
	vector<double> sigma( sigma_start ); 
	vector<double> sigma_last( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_start_i( dim ); 

	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used

//	double temp, max_diff = 1.0;
    double mean_diff = 1.0;
	while( mean_diff > threshold_c )
	{
		memcpy( &sigma_last[0], &sigma[0], sizeof( double ) * pxp );
		
		for( i = 0; i < dim; i++ )
		{
			// Count size of note
			size_node = 0;
			for( j = 0; j < dim; j++ ) size_node += G[ j * dim + i ];

			if( size_node > 0 )
			{
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					if( G[ j * dim + i ] )
					{
						sigma_start_N_i[ l ] = sigma_start[ i * dim + j ]; 
						N_i[ l++ ]           = j;
					}else
						beta_star[ j ] = 0.0; 
				}
				// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
				
				sub_matrix( &sigma[0], &sigma_N_i[0], &N_i[0], &size_node, &dim );
					
				// A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
				F77_NAME(dposv)( &uplo, &size_node, &one, &sigma_N_i[0], &size_node, &sigma_start_N_i[0], &size_node, &info FCONE );

				for( j = 0; j < size_node; j++ ) beta_star[ N_i[ j ] ] = sigma_start_N_i[ j ];
				
				F77_NAME(dgemm)( &transN, &transN, &dim, &one, &dim, &alpha, &sigma[0], &dim, &beta_star[0], &dim, &beta, &sigma_start_i[0], &dim FCONE FCONE );
				
				for( j = 0; j < i; j++ )
				{
					sigma[ j * dim + i ] = sigma_start_i[ j ];
					sigma[ i * dim + j ] = sigma_start_i[ j ];
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					sigma[ j * dim + i ] = sigma_start_i[ j ];
					sigma[ i * dim + j ] = sigma_start_i[ j ];
				}
			}else{
				for( j = 0; j < i; j++ )
				{
					sigma[ j * dim + i ] = 0.0;
					sigma[ i * dim + j ] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					sigma[ j * dim + i ] = 0.0;
					sigma[ i * dim + j ] = 0.0;
				}
			} 
		}

		mean_diff = fabs( static_cast<double>( sigma[ 0 ] - sigma_last[ 0 ] ) );
		for( i = 1; i < pxp; i++ )
		    mean_diff += fabs( static_cast<double>( sigma[ i ] - sigma_last[ i ] ) );
		mean_diff /= pxp;

//		max_diff = fabs( static_cast<double>( sigma[ 0 ] - sigma_last[ 0 ] ) );
//		for( i = 1; i < pxp; i++ )
//		{
//			temp = fabs( static_cast<double>( sigma[ i ] - sigma_last[ i ] ) );
//			if( temp > max_diff ) max_diff = temp; 
//		}		
	}

	inverse( &sigma[0], K, &dim );
}
     
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// rgwish ONLY for inside of MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rgwish_sigma( int G[], int size_node[], double Ts[], double K[], double sigma[], int *b_star, 
					int *p, double *threshold,
					double sigma_start[], double inv_C[], double beta_star[], double sigma_i[], 
					vector<double> &sigma_start_N_i, vector<double> &sigma_N_i, vector<int> &N_i )
{
	int i, i1, j, ij, ip, l, size_node_i, info, one = 1, dim = *p, pxp = dim * dim, dim1 = dim + 1, bKdim = *b_star + dim - 1;	
	
	double alpha = 1.0, beta  = 0.0;	
	
	char transT  = 'T', transN = 'N', side = 'R', upper = 'U';																	
	
	// - - STEP 1: sampling from wishart distributions  - - - - - - - - - - - - - - - - - - - - - -|
	// - -  Sample values in Psi matrix
	//GetRNGstate();
	#pragma omp parallel for
	for( i = 0; i < dim; i++ )
		sigma_start[ i * dim1 ] = sqrt( Rf_rgamma( ( bKdim - i ) * 0.5, 2.0 ) ); // i * dim1 = i * dim + i
		//sigma_start[i * dim1] = sqrt( rchisq( bKdim - i ) ); // i * dim1 = i * dim + i

	#pragma omp parallel for
	for( j = 1; j < dim; j++ )
		for( int i = 0; i < j; i++ )
		{
			sigma_start[ j * dim + i ] = norm_rand();
			sigma_start[ i * dim + j ] = 0.0;
		}
	//PutRNGstate();
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	
	// C = psi %*% Ts   I used psi = psi %*% Ts   Now is  sigma_start = sigma_start %*% Ts
	F77_NAME(dtrmm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, Ts, &dim, &sigma_start[0], &dim FCONE FCONE FCONE FCONE );

	side = 'L';
	// creating an identity matrix
	#pragma omp parallel for
	for( i = 0; i < dim; i++ )
		for( int j = 0; j < dim; j++ )
			inv_C[ j * dim + i ] = ( i == j );	
	
	// op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
	F77_NAME(dtrsm)( &side, &upper, &transN, &transN, &dim, &dim, &alpha, &sigma_start[0], &dim, &inv_C[0], &dim FCONE FCONE FCONE FCONE );
 
	// sigma_start = inv_C %*% t( inv_C )  																				
	F77_NAME(dgemm)( &transN, &transT, &dim, &dim, &dim, &alpha, &inv_C[0], &dim, &inv_C[0], &dim, &beta, &sigma_start[0], &dim FCONE FCONE );

	memcpy( sigma, &sigma_start[0], sizeof( double ) * pxp ); 
	
//	double temp, max_diff = 1.0, threshold_c = *threshold;
	double mean_diff = 1.0, threshold_c = *threshold;
    int counter = 0;
	while( ( mean_diff > threshold_c ) and ( counter < 5000 ) )
	{
		counter++;
	    mean_diff = 0.0;
		
		for( i = 0; i < dim; i++ )
		{
			ip = i * dim;

			size_node_i = size_node[ i ];
			if( size_node_i > 0 )
			{
				l = 0;
				for( j = 0; j < dim; j++ )
				{
					ij = ip + j;
					if( G[ ij ] )
					{
						sigma_start_N_i[ l ] = sigma_start[ ij ]; 
						N_i[ l++ ]           = j;
					}
					else
						beta_star[ j ] = 0.0; 
				}
				// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
				
				sub_matrix_upper( sigma, &sigma_N_i[0], &N_i[0], &size_node_i, &dim );
					
				// A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
				F77_NAME(dposv)( &upper, &size_node_i, &one, &sigma_N_i[0], &size_node_i, &sigma_start_N_i[0], &size_node_i, &info FCONE );

				for( j = 0; j < size_node_i; j++ ) beta_star[ N_i[ j ] ] = sigma_start_N_i[ j ];
	
				// sigma_i = sigma %*% beta_star
				F77_NAME(dsymv)( &side, &dim, &alpha, sigma, &dim, &beta_star[0], &one, &beta, &sigma_i[0], &one FCONE );
				
				memcpy( sigma + ip, sigma_i, sizeof( double ) * i );	
				
				for( j = 0; j < i; j++ )
				{
					ij         = j * dim + i;
					mean_diff += fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//					temp      = fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//					max_diff  = ( temp > max_diff ) ? temp : max_diff; 					

					sigma[ ij ] = sigma_i[ j ];
				}
				
				i1 = i + 1;
				memcpy( sigma + ip + i1, sigma_i + i1, sizeof( double ) * ( dim - i1 ) );	

				for( j = i1; j < dim; j++ )
				{
					ij         = j * dim + i;
					mean_diff += fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//					temp      = fabs( static_cast<double>( sigma[ ij ] - sigma_i[ j ] ) );
//					max_diff  = ( temp > max_diff ) ? temp : max_diff; 					

					sigma[ ij ] = sigma_i[ j ];
				}
			}else{					
				for( j = 0; j < i; j++ )
				{
					ij         = j * dim + i;
				    mean_diff += fabs( static_cast<double>( sigma[ ij ] ) );
//				    temp     = fabs( static_cast<double>( sigma[ ij ] ) );
//					max_diff = ( temp > max_diff ) ? temp : max_diff; 					

					sigma[ ij ]     = 0.0;
					sigma[ ip + j ] = 0.0;
				}
				
				for( j = i + 1; j < dim; j++ )
				{
					ij         = j * dim + i;
					mean_diff += fabs( static_cast<double>( sigma[ ij ] ) );
//					temp     = fabs( static_cast<double>( sigma[ ij ] ) );
//					max_diff = ( temp > max_diff ) ? temp : max_diff; 					

					sigma[ ij     ] = 0.0;
					sigma[ ip + j ] = 0.0;				
				}
			} 
		}
		
		mean_diff /= pxp;
	}
	
	memcpy( &sigma_start[0], sigma, sizeof( double ) * pxp );	 	
	
	inverse( &sigma_start[0], K, &dim );
	// creating an identity matrix
	//#pragma omp parallel for
	//for( i = 0; i < dim; i++ )
	//	for( int j = 0; j < dim; j++ )
	//		K[ j * dim + i ] = ( i == j );
	
	// LAPACK function: computes solution to A * X = B, where A is symmetric positive definite matrix
	//F77_NAME(dposv)( &upper, &dim, &dim, &sigma_start[0], &dim, K, &dim, &info );
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Part of function "gnorm"
// which is for calculating Normalizing constant of G-Wishart distribution 
// based on Monto Carlo algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_exp_mc( int G[], int nu[], int *b, double H[], int *check_H, int *mc, int *p, double f_T[] )
{
	int iter, i, j, ij, h, r, mc_iter = *mc, dim = *p, pxp = dim * dim, b_c = *b;
	
	double sumPsi, sumPsiH, sumPsiHi, sumPsiHj;
	double max_numeric_limits_ld = numeric_limits<double>::max() / 1000;
	double min_numeric_limits_ld = numeric_limits<double>::min() * 1000;
	
	vector<double> psi( pxp, 0.0 );      

	GetRNGstate();
	if( *check_H == 1 )
	{ 
		for( iter = 0; iter < mc_iter; iter++ ) 
		{
			for( i = 0; i < dim; i++ )
				psi[ i * dim + i ] = sqrt( Rf_rgamma( ( b_c + nu[ i ] ) / 2.0, 2.0 ) );
				//psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );

			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					//if( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 ); else psi[ij] = 0.0;
					psi[ ij ] = ( G[ ij ] == 1 ) ? norm_rand() : 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if( G[ ij ] == 0 )
					{
						psi[ ij ] = 0.0; // it's not necessary
						if( i > 0 )  
						{
							sumPsi = 0.0;
							//sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] )
							// for( h = 0; h < ( i - 1 ); h++ )
							for( h = 0; h < i; h++ )
							{
								if( sumPsi == R_PosInf ) sumPsi = max_numeric_limits_ld;	
								if( sumPsi == R_NegInf ) sumPsi = min_numeric_limits_ld;	
								sumPsi += ( psi[ i * dim + h ] * psi[ j * dim + h ] );
							}
							
							//psi[i, j] <- - sum( psi[ 1 : ( i - 1 ), i ] * psi[ 1 : ( i - 1 ), j ] ) / psi[i, i]
							psi[ ij ] = - sumPsi / psi[ i * dim + i ];
						}
						
						if( psi[ ij ] == R_PosInf ) psi[ ij ] = max_numeric_limits_ld;	
						if( psi[ ij ] == R_NegInf ) psi[ ij ] = min_numeric_limits_ld;	
												
						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[ iter ] += ( psi[ ij ] * psi[ ij ] ); 
					}
				}
		
			// checking Inf values
			if( f_T[ iter ] == R_PosInf ) f_T[ iter ] = max_numeric_limits_ld;			
		} 
	}else{
		for( iter = 0; iter < mc_iter; iter++ ) 
		{
			for( i = 0; i < dim; i++ )
				psi[ i * dim + i ] = sqrt( Rf_rgamma( ( b_c + nu[ i ] ) / 2.0, 2.0 ) );
				//psi[i * dim + i] = sqrt( rchisq( b_c + nu[i] ) );

			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					//if( G[ij] == 1 ) psi[ij] = rnorm( 0, 1 ); elsepsi[ij] = 0.0;
					psi[ ij ] = ( G[ ij ] == 1 ) ? norm_rand() : 0.0;
				}
			
			for( i = 0; i < dim - 1; i++ ) 
				for( j = i + 1; j < dim; j++ ) 
				{
					ij = j * dim + i;
					
					if( G[ ij ] == 0 )
					{
						//psi[i, j] = - sum( psi[ i, i : ( j - 1 ) ] * H[ i : ( j - 1 ), j ] )
						sumPsiH = 0.0;
						for( h = i; h < j; h++ )
						{
							if( sumPsiH == R_PosInf ) sumPsiH = max_numeric_limits_ld;	
							if( sumPsiH == R_NegInf ) sumPsiH = min_numeric_limits_ld;	
							sumPsiH += ( psi[ h * dim + i ] * H[ j * dim + h ] ); 
						}
						psi[ ij ] = - sumPsiH;
						
						if( i > 0 )  //if( i > 1 )
							for( r = 0; r < i; r++ ) //for( r in 1 : ( i - 1 ) )
							{
								//sum( psi[ r, r : i ] * H[ r : i, i ] )
								sumPsiHi = 0.0;
								for( h = r; h < i + 1; h++  )
								{
									if( sumPsiHi == R_PosInf ) sumPsiHi = max_numeric_limits_ld;	
									if( sumPsiHi == R_NegInf ) sumPsiHi = min_numeric_limits_ld;	
									sumPsiHi += ( psi[ h * dim + r ] * H[ i * dim + h ] );	
								}
									
								//sum( psi[ r, r : j ] * H[ r : j, j ] ) )
								sumPsiHj = 0.0;
								for( h = r; h < j + 1; h++  )
									sumPsiHj += ( psi[ h * dim + r ] * H[ j * dim + h ] );
								
								//psi[i, j] <- psi[i, j] - ( ( sum( psi[ r, r : i ] * H[ r : i, i ] ) ) * ( sum( psi[ r, r : j ] * H[ r : j, j ] ) ) ) / ( psi[i, i] )
								psi[ ij ] -= ( sumPsiHi * sumPsiHj ) / psi[ i * dim + i ];
							}

						if( psi[ ij ] == R_PosInf ) psi[ ij ] = max_numeric_limits_ld;	
						if( psi[ ij ] == R_NegInf ) psi[ ij ] = min_numeric_limits_ld;	

						//f_T[k] <- f_T[k] + psi[i, j] ^ 2
						f_T[ iter ] += ( psi[ ij ] * psi[ ij ] ); 
					}
				}
		
			// checking Inf values
			if( f_T[ iter ] == R_PosInf ) f_T[ iter ] = max_numeric_limits_ld;			
		}
	}
	PutRNGstate();	
} 
      
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
     
