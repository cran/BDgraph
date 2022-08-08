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

#include "matrix.h"
#include "rgwish.h"
#include "copula.h"

extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing alpha (probability of acceptance) in RJ-MCMC algorithm
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_alpha_rjmcmc( double *log_alpha_ij, double log_ratio_g_prior[], int *selected_edge_i, 
                    int *selected_edge_j, int G[], double Ds[], 
                    double sigma[], double sigma21[], double sigma22[], double sigmaj12[], double sigmaj22[],    
                    double K[], double K21[], double K121[], double Kj12[], 
                    double K12xK22_inv[], double Kj12xK22_inv[], double sigma11_inv[], 
                    double sigma21xsigma11_inv[], int *b, int *p )
{
	int one = 1, two = 2, dim = *p, dim1 = dim + 1, p1 = dim - 1, p2 = dim - 2;

	double alpha = 1.0, beta = 0.0, alpha1 = -1.0, beta1 = 1.0;
	char transT = 'T', transN = 'N', sideL = 'L';																	

	int ij      = *selected_edge_j * dim + *selected_edge_i;
	int jj      = *selected_edge_j * dim1;
	double Dsij = Ds[ ij ];
	double Dsjj = Ds[ jj ];
	
	sub_matrices1( &sigma[0], &sigmaj12[0], &sigmaj22[0], selected_edge_j, &dim );

	// sigma[-j,-j] - ( sigma[-j, j] %*% sigma[j, -j] ) / sigma[j,j]
	// Kj22_inv <- sigmaj22 = sigmaj22 - sigmaj12 * sigmaj12 / sigmaj11
	double sigmajj_inv = - 1.0 / sigma[jj];
	F77_NAME(dsyr)( &sideL, &p1, &sigmajj_inv, &sigmaj12[0], &one, &sigmaj22[0], &p1 FCONE );
		
	// For (i,j) = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
	sub_row_mins( K, &Kj12[0], selected_edge_j, &dim );   // Kj12 = K[j, -j]  
	Kj12[ *selected_edge_i ] = 0.0;                         // Kj12[1,i] = 0

	// Kj12xK22_inv = Kj12 %*% Kj22_inv here sigmaj22 instead of Kj22_inv
	F77_NAME(dsymv)( &sideL, &p1, &alpha, &sigmaj22[0], &p1, &Kj12[0], &one, &beta, &Kj12xK22_inv[0], &one FCONE );
	
	// K022 = Kj12xK22_inv %*% t(Kj12)
	double K022 = F77_NAME(ddot)( &p1, &Kj12xK22_inv[0], &one, &Kj12[0], &one );			

	// For (i,j) = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	sub_cols_mins( K, &K21[0], selected_edge_i, selected_edge_j, &dim );  // K21 = K[-e, e]  
	
	sub_matrices_inv( &sigma[0], &sigma11_inv[0], &sigma21[0], &sigma22[0], selected_edge_i, selected_edge_j, &dim );

	// sigma21xsigma11_inv = sigma21 %*% sigma11_inv
	F77_NAME(dgemm)( &transN, &transN, &p2, &two, &two, &alpha, &sigma21[0], &p2, &sigma11_inv[0], &two, &beta, &sigma21xsigma11_inv[0], &p2 FCONE FCONE );

	// sigma22 = sigma22 - sigma21xsigma11_inv %*% t( sigma21 )
	F77_NAME(dgemm)( &transN, &transT, &p2, &p2, &two, &alpha1, &sigma21xsigma11_inv[0], &p2, &sigma21[0], &p2, &beta1, &sigma22[0], &p2 FCONE FCONE );

	// K12xK22_inv = t( K21 ) %*% K22_inv  here sigam12 = K22_inv
	F77_NAME(dgemm)( &transT, &transN, &two, &p2, &p2, &alpha, &K21[0], &p2, &sigma22[0], &p2, &beta, &K12xK22_inv[0], &two FCONE FCONE );  
	
	// K121 = K12xK22_inv %*% K21													
	F77_NAME(dgemm)( &transN, &transN, &two, &two, &p2, &alpha, &K12xK22_inv[0], &two, &K21[0], &p2, &beta, &K121[0], &two FCONE FCONE );		
	// Finished (i,j) = 1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

	double a11      = K[*selected_edge_i * dim1] - K121[0];	
	double sum_diag = Dsjj * ( K022 - K121[3] ) - Dsij * ( K121[1] + K121[2] );

	// nu_star = b + sum( Gf[,i] * Gf[,j] )
	int nu_star = *b;
	for( int k = 0; k < dim; k++ )
		nu_star += G[*selected_edge_i * dim + k] * G[*selected_edge_j * dim + k];
	nu_star = 0.5 * nu_star;

	// *log_alpha_ij = 0.5 * ( log( static_cast<double>( 2.0 ) ) + log( static_cast<double>( Dsjj ) ) - log( static_cast<double>( a11 ) ) ) + 
	*log_alpha_ij = - log_ratio_g_prior[ ij ] + 0.5 * log( 2.0 * Dsjj / a11 ) + 
			  lgammafn( nu_star + 0.5 ) - lgammafn( nu_star ) - 0.5 * ( Dsij * Dsij * a11 / Dsjj  + sum_diag );

	if( G[ ij ] == 0 ) *log_alpha_ij = - *log_alpha_ij;	
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian Graphical models  
// for D = I_p 
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_rjmcmc_ma( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double K[], 
                    int *p, double *threshold, double K_hat[], int p_links[], 
                    int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin;
	int selected_edge, selected_edge_i, selected_edge_j, ip, i, j, ij, counter;
	int dim = *p, pxp = dim * dim, p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2;
	int qp = dim * ( dim - 1 ) / 2;
	double log_alpha_ij;

	// - - allocation for log_alpha_ij 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// - -  for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// - - - - - - - - - - - - - - 

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of selected edge 
	vector<int> index_row( qp );
	vector<int> index_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = g_prior[ j * dim + i ];
	        if( ( ij != 0.0 ) or ( ij != 1.0 ) )
	        {
	            index_row[ counter ] = i;
	            index_col[ counter ] = j;
	            counter++;
	        }
	    }
	int sub_qp = counter;

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + 1 ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	  				
// - - - STEP 1: selecting edge and calculating alpha - - - - - - - - - - - - - - - - - - - - - - -|		
		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		selected_edge   = static_cast<int>( unif_rand() * sub_qp );
		selected_edge_i = index_row[ selected_edge ];
		selected_edge_j = index_col[ selected_edge ];

// - - - STEP 1: calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		

		log_alpha_rjmcmc( &log_alpha_ij, &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, G, Ds, 
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );

// - - - End of calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |		
		  		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < log_alpha_ij )
		{
			ij    = selected_edge_j * dim + selected_edge_i;
			G[ ij ] = 1 - G[ ij ];
			G[ selected_edge_i * dim + selected_edge_j ] = G[ ij ];

			if( G[ ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{ 
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}
		}

// - - - STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
			for( i = 0; i < pxp ; i++ )
			{
				K_hat[ i ]   += K[ i ];
				p_links[ i ] += G[ i ];
			}	
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
	}  
	PutRNGstate();
// - - - End of main loop for reversible jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - | 
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian Graphical models  
// for D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_rjmcmc_map( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double K[], 
                    int *p, double *threshold, 
                    int all_graphs[], double all_weights[], double K_hat[], 
                    char *sample_graphs[], double graph_weights[], int *size_sample_g,
                    int *b, int *b_star, double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, count_all_g = 0;
	int selected_edge, selected_edge_i, selected_edge_j, size_sample_graph = *size_sample_g;
	int ip, i, j, ij, counter, dim = *p, pxp = dim * dim, p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2;
	int qp = dim * ( dim - 1 ) / 2;

	double log_alpha_ij;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );

	bool this_one;

	// - - - allocation for log_alpha_ij 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// - - for rgwish_sigma 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// - - - - - - - - - - - - - - 
	vector<char> char_g( qp );               // char string_g[pp];

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	vector<int> index_row( qp );
	vector<int> index_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = g_prior[ j * dim + i ];
	        if( ( ij != 0.0 ) or ( ij != 1.0 ) )
	        {
	            index_row[ counter ] = i;
	            index_col[ counter ] = j;
	            counter++;
	        }
	    }
	int sub_qp = counter;

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - - Main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + 1 ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	  				
		// STEP 1: selecting edge and calculating alpha
		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		selected_edge   = static_cast<int>( unif_rand() * sub_qp );
		selected_edge_i = index_row[ selected_edge ];
		selected_edge_j = index_col[ selected_edge ];

// - - - STEP 1: calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		

		log_alpha_rjmcmc( &log_alpha_ij, &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, G, Ds, 
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );

// - - - End of calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |		
		  		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < log_alpha_ij )
		{
			ij    = selected_edge_j * dim + selected_edge_i;
			G[ ij ] = 1 - G[ ij ];
			G[ selected_edge_i * dim + selected_edge_j ] = G[ ij ];

			if( G[ ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{ 
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}
		}

// - - - STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			counter = 0;	
			for( j = 1; j < dim; j++ )
				for( i = 0; i < j; i++ )
					char_g[ counter++ ] = G[ j * dim + i ] + '0'; 
		
			for( i = 0; i < pxp ; i++ ) K_hat[ i ] += K[ i ];	

			string_g = string( char_g.begin(), char_g.end() );	
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[ i ] == string_g )
				{
					graph_weights[ i ]++;           // += all_weights[count_all_g];
					all_graphs[ count_all_g ] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[ size_sample_graph ] = string_g;
				graph_weights[ size_sample_graph ]   = all_weights[ count_all_g ];
				all_graphs[ count_all_g ]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
	}  
	PutRNGstate();
// - - - End of main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
               
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian copula Graphical models  
// for D = I_p 
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_rjmcmc_ma( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double K[], 
                    int *p, double *threshold, 
                    double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                    double K_hat[], int p_links[], 
                    int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin;
	int selected_edge, counter, selected_edge_i, selected_edge_j;
	int ip, i, j, ij, dim = *p, pxp = dim * dim, p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2;
	int qp = dim * ( dim - 1 ) / 2;
	
	double log_alpha_ij;

	// - - allocation for log_alpha_ij 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// - - for rgwish_sigma  - - - - - - - - - 
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// - - for copula  - - - - - - - - - - - - 
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// - - - - - - - - - - - - - - - - - - - -

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of selected edge 
	vector<int> index_row( qp );
	vector<int> index_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = g_prior[ j * dim + i ];
	        if( ( ij != 0.0 ) or ( ij != 1.0 ) )
	        {
	            index_row[ counter ] = i;
	            index_col[ counter ] = j;
	            counter++;
	        }
	    }
	int sub_qp = counter;
	
	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - - Main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + 1 ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	  		
// - - - STEP 1: copula - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
		
		get_Ds( K, Z, R, not_continuous, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
// - - - STEP 2: calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		selected_edge   = static_cast<int>( unif_rand() * sub_qp );
		selected_edge_i = index_row[ selected_edge ];
		selected_edge_j = index_col[ selected_edge ];
		
// - - - STEP 2: calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		

		log_alpha_rjmcmc( &log_alpha_ij, &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, G, Ds, 
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );

// - - - End of calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |		
		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < log_alpha_ij )
		{
			ij    = selected_edge_j * dim + selected_edge_i;
			G[ ij ] = 1 - G[ ij ];
			G[ selected_edge_i * dim + selected_edge_j ] = G[ ij ];

			if( G[ ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{ 
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}
		}

// - - - STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
			for( i = 0; i < pxp ; i++ )
			{
				K_hat[ i ]   += K[ i ];
				p_links[ i ] += G[ i ];
			}	
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
	}  
	PutRNGstate();
// - - - End of main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - | 
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian copula Graphical models  
// for D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_rjmcmc_map( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double K[], 
                    int *p, double *threshold, 
                    double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                    int all_graphs[], double all_weights[], double K_hat[], 
                    char *sample_graphs[], double graph_weights[], int *size_sample_g,
                    int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, count_all_g = 0;
	int selected_edge, counter, selected_edge_i, selected_edge_j, size_sample_graph = *size_sample_g;
	int ip, i, j, ij, dim = *p, pxp = dim * dim, p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2;
	int qp = dim * ( dim - 1 ) / 2;
	bool this_one;
	double log_alpha_ij;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	// - - - allocation for log_alpha_ij 
	vector<double> K121( 4 ); 
	vector<double> Kj12( p1 );              // K[j, -j]
	vector<double> sigmaj12( p1 );          // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );     // sigma[-j, -j]
	vector<double> Kj12xK22_inv( p1 ); 

	vector<double> K21( p2x2 );                          // K[-e, e]
	vector<double> sigma21( p2x2 );                      // sigma[-e, e]
	vector<double> sigma22( ( dim - 2 ) * ( dim - 2 ) ); // sigma[-e, -e]
	vector<double> sigma11_inv( 4 );                     // inv( sigma[e, e] )
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> K12xK22_inv( p2x2 );   
	// - - for rgwish_sigma - - - - - - - - -
	vector<double> sigma_start( pxp ); 
	vector<double> inv_C( pxp ); 
	vector<double> beta_star( dim ); 
	vector<double> sigma_i( dim ); 
	vector<double> sigma_start_N_i( dim );   // For dynamic memory used
	vector<double> sigma_N_i( pxp );         // For dynamic memory used
	vector<int> N_i( dim );                  // For dynamic memory used
	// - - for copula - - - - - - - - - - - - 
	vector<double> S( pxp ); 
	vector<double> inv_Ds( pxp ); 
	vector<double> copy_Ds( pxp ); 
	// - - - - - - - - - - -- - - - - - - - -
	vector<char> char_g( qp );              // char string_g[pp];
	
	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of selected edge 
	vector<int> index_row( qp );
	vector<int> index_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = g_prior[ j * dim + i ];
	        if( ( ij != 0.0 ) or ( ij != 1.0 ) )
	        {
	            index_row[ counter ] = i;
	            index_col[ counter ] = j;
	            counter++;
	        }
	    }
	int sub_qp = counter;

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + 1 ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	  		
// - - - STEP 1: copula - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
		
		get_Ds( K, Z, R, not_continuous, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
// - - - STEP 2: calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		selected_edge   = static_cast<int>( unif_rand() * sub_qp );
		selected_edge_i = index_row[ selected_edge ];
		selected_edge_j = index_col[ selected_edge ];

// - - - STEP 2: calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		

		log_alpha_rjmcmc( &log_alpha_ij, &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, G, Ds, 
                      &sigma[0], &sigma21[0], &sigma22[0], &sigmaj12[0], &sigmaj22[0],    
                      &K[0], &K21[0], &K121[0], &Kj12[0], 
                      &K12xK22_inv[0], &Kj12xK22_inv[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],  
                      b, &dim );

// - - - End of calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |		
		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < log_alpha_ij )
		{
			ij    = selected_edge_j * dim + selected_edge_i;
			G[ ij ] = 1 - G[ ij ];
			G[ selected_edge_i * dim + selected_edge_j ] = G[ ij ];

			if( G[ ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{ 
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}
		}

// - - - STEP 2: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			counter = 0;	
			for( j = 1; j < dim; j++ )
				for( i = 0; i < j; i++ )
					char_g[ counter++ ] = G[ j * dim + i ] + '0'; 
		
			for( i = 0; i < pxp ; i++ ) K_hat[ i ] += K[ i ];	

			string_g = string( char_g.begin(), char_g.end() );	
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[ i ] == string_g )
				{
					graph_weights[ i ]++;     // += all_weights[count_all_g];
					all_graphs[ count_all_g ] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[ size_sample_graph ] = string_g;
				graph_weights[ size_sample_graph ]   = all_weights[ count_all_g ];
				all_graphs[ count_all_g ]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
	}  
	PutRNGstate();
// - - End of main loop for Reversible Jump MCMC - - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy(sample_graphs[ i ], qp, 0);
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
    
} // End of exturn "C"
