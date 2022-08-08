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
// birth-death MCMC for Gaussian copula Graphical models  
// Based on Double Metropolis-Hastings
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_DMH_bdmcmc_ma( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double Ti[], 
                        double K[], int *p, double *threshold, 
                        double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                        double K_hat[], double p_links[],
                        int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, b1 = *b;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, counter;
	int ip, i, j, ij, dim = *p, pxp = dim * dim, one = 1;
	int qp = dim * ( dim - 1 ) / 2;
	double sum_weights = 0.0, sum_rates, weight_C;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	vector<double> K121( 4 ); 															

	// - -  for rgwish_sigma - - - - - - - - - 
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
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of rates 
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
	vector<double> rates( sub_qp );

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
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
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		rates_bdmcmc_dmh_parallel( &rates[0], &log_ratio_g_prior[0], &G[0], &index_row[0], &index_col[0], &sub_qp, &Ds[0], &D[0], &sigma[0], &K[0], &sigma_dmh[0], &K_dmh[0], b, &dim );
		
// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			for( i = 0; i < pxp ; i++ )
				if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
			
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[ selected_edge_ij ] = 1 - G[ selected_edge_ij ];
		G[ selected_edge_i * dim + selected_edge_j ] = G[ selected_edge_ij ];

		if( G[ selected_edge_ij ] )
		{ 
			++size_node[ selected_edge_i ]; 
			++size_node[ selected_edge_j ]; 
		}
		else
		{ 
			--size_node[ selected_edge_i ]; 
			--size_node[ selected_edge_j ]; 
		}

// - - - STEP 3: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < pxp; i++ )
	{	
		p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
		K_hat[ i ]   = K_hat_Cpp[ i ] / sum_weights;
	}
}
   
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian copula Graphical models  
// Based on Double Metropolis-Hastings
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_DMH_bdmcmc_map( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double Ti[], 
                        double K[], int *p, double *threshold, 
                        double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                        int all_graphs[], double all_weights[], double K_hat[], 
                        char *sample_graphs[], double graph_weights[], int *size_sample_g,
                        int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, b1 = *b, count_all_g = 0;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, counter, size_sample_graph = *size_sample_g;
	int ip, i, j, ij, dim = *p, pxp = dim * dim, one = 1;
	int qp = dim * ( dim - 1 ) / 2;
	
	double sum_weights = 0.0, sum_rates, weight_C;
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<char> char_g( qp );              // char string_g[pp];
	vector<double> K121( 4 ); 																

	// - -  for rgwish_sigma  - - - - - - - - -
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
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	// - - - - - - - - - - - - - - - - - - - - 

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of rates 
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
	vector<double> rates( sub_qp );

	vector<double>log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
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
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		rates_bdmcmc_dmh_parallel( &rates[0], &log_ratio_g_prior[0], &G[0], &index_row[0], &index_col[0], &sub_qp, &Ds[0], &D[0], &sigma[0], &K[0], &sigma_dmh[0], &K_dmh[0], b, &dim );

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[ count_all_g ] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[ i ] == string_g )
				{
					graph_weights[ i ] += all_weights[ count_all_g ];
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
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[ selected_edge_ij ] = 1 - G[ selected_edge_ij ];
		G[ selected_edge_i * dim + selected_edge_j ] = G[ selected_edge_ij ];

		if( G[ selected_edge_ij ] )
		{ 
			++size_node[ selected_edge_i ]; 
			++size_node[ selected_edge_j ]; 
		}else{ 
			--size_node[ selected_edge_i ]; 
			--size_node[ selected_edge_j ]; 
		}

// - - - STEP 3: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;

	for( i = 0; i < pxp; i++ ) 
		K_hat[ i ] /= sum_weights;
}
           
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Multiple birth-death MCMC for Gaussian copula Graphical models  
// Based on Double Metropolis-Hastings
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_DMH_bdmcmc_ma_multi_update( int *iter, int *burnin, int G[], double g_prior[], double Ts[], 
                        double Ti[], double K[], int *p, double *threshold, 
                        double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                        double K_hat[], double p_links[],
                        int *b, int *b_star, double D[], double Ds[], int *multi_update, int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, b1 = *b, multi_update_C = *multi_update;
	int selected_edge_i, selected_edge_j, selected_edge_ij, counter;
	int ip, i, j, ij, dim = *p, pxp = dim * dim, one = 1;
	int qp = dim * ( dim - 1 ) / 2;
	double sum_weights = 0.0, sum_rates, weight_C; 	

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	vector<double> K121( 4 ); 																

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
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of rates 
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
	vector<double> rates( sub_qp );

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + size_index ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + size_index ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}

// - - - STEP 1: copula - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
		
		get_Ds( K, Z, R, not_continuous, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		rates_bdmcmc_dmh_parallel( &rates[0], &log_ratio_g_prior[0], &G[0], &index_row[0], &index_col[0], &sub_qp, &Ds[0], &D[0], &sigma[0], &K[0], &sigma_dmh[0], &K_dmh[0], b, &dim );

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			for( i = 0; i < pxp ; i++ )
				if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
			
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &sub_qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_row[ index_selected_edges[ i ] ];
			selected_edge_j = index_col[ index_selected_edges[ i ] ];
			
			selected_edge_ij      = selected_edge_j * dim + selected_edge_i;
			G[ selected_edge_ij ] = 1 - G[ selected_edge_ij ];
			G[ selected_edge_i * dim + selected_edge_j ] = G[ selected_edge_ij ];
		
			if( G[ selected_edge_ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{ 
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}		
		}
		
// - - - STEP 3: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < pxp; i++ )
	{	
		p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
		K_hat[ i ]   = K_hat_Cpp[ i ] / sum_weights;
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Multiple birth-death MCMC for Gaussian copula Graphical models  
// Based on Double Metropolis-Hastings
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_DMH_bdmcmc_map_multi_update( int *iter, int *burnin, int G[], double g_prior[], double Ts[], 
                    double Ti[], double K[], int *p, double *threshold, 
                    double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                    int all_graphs[], double all_weights[], double K_hat[], 
                    char *sample_graphs[], double graph_weights[], int *size_sample_g, int *counter_all_g,
                    int *b, int *b_star, double D[], double Ds[], int *multi_update, int *print )
{
	int print_c = *print, multi_update_C = *multi_update, iteration = *iter, burn_in = *burnin, b1 = *b;
	int count_all_g = *counter_all_g;
	int selected_edge_i, selected_edge_j, selected_edge_ij, counter, size_sample_graph = *size_sample_g;
	int ip, i, j, ij, dim = *p, pxp = dim * dim, one = 1;
	int qp = dim * ( dim - 1 ) / 2;
	double sum_weights = 0.0, sum_rates, weight_C; 
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<char> char_g( qp );              // char string_g[pp];
	vector<double> K121( 4 ); 																

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
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings
	// - - - - - - - - - - - - - - - - - - - - 

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of rates 
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
	vector<double> rates( sub_qp );

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + size_index ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + size_index ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}

// - - - STEP 1: copula - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|		
		
		get_Ds( K, Z, R, not_continuous, D, Ds, &S[0], gcgm, n, &dim );
		get_Ts( Ds, Ts, &inv_Ds[0], &copy_Ds[0], &dim );
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		rates_bdmcmc_dmh_parallel( &rates[0], &log_ratio_g_prior[0], &G[0], &index_row[0], &index_col[0], &sub_qp, &Ds[0], &D[0], &sigma[0], &K[0], &sigma_dmh[0], &K_dmh[0], b, &dim );

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			//for( i = 0; i < pxp; i++ ) K_hat[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat[0], &one );			

			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[ count_all_g ] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[ i ] == string_g )
				{
					graph_weights[ i ] += all_weights[ count_all_g ];
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
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &sub_qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_row[ index_selected_edges[ i ] ];
			selected_edge_j = index_col[ index_selected_edges[ i ] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[ selected_edge_ij ] = 1 - G[ selected_edge_ij ];
			G[ selected_edge_i * dim + selected_edge_j ] = G[ selected_edge_ij ];
		
			if( G[ selected_edge_ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{ 
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}		
		}

// - - - STEP 3: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy(sample_graphs[ i ], qp, 0);
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	*counter_all_g = count_all_g;

	for( i = 0; i < pxp; i++ ) K_hat[ i ] /= sum_weights;
}
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian copula Graphical models  
// Based on Double Metropolis-Hastings
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_DMH_rjmcmc_ma( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double Ti[], 
                        double K[], int *p, double *threshold, 
                        double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                        double K_hat[], double p_links[],
                        int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, b1 = *b;
	int index_selected_edge, selected_edge_i, selected_edge_j, counter;
	int ip, i, j, ij, jj, dim = *p, pxp = dim * dim;
	int qp = dim * ( dim - 1 ) / 2;
	int p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2;
	double Dsijj, Dsjj, Dsij, logH_ij, logI_p, Dij, Djj, Dijj, alpha_ij;   


	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	vector<double> K121( 4 ); 															
	
	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );      // sigma[-j, -j]
	vector<double> Kj22_inv( p1 * p1 ); 
	vector<double> Kj12xK22_inv( p1 ); 
	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2 * p2 );       // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2 * p2 ); 
	vector<double> K22_inv( p2 * p2 ); 
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
	// - - - - - - - - - - - - - - - - - - - - 
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of rates 
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
	//vector<double> rates( sub_qp );

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
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
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		index_selected_edge   = static_cast<int>( unif_rand() * sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

		// - - - -  Calculating alpha - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
		ij    = selected_edge_j * dim + selected_edge_i;
		jj    = selected_edge_j * dim + selected_edge_j;
		Dsij  = Ds[ ij ];
		Dsjj  = Ds[ jj ];
		Dsijj = - Dsij * Dsij / Dsjj;
		Dij   = D[ ij ];
		Djj   = D[ jj ];
		Dijj  = - Dij * Dij / Djj;

		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		log_H_ij( &K[0], &sigma[0], &logH_ij, &selected_edge_i, &selected_edge_j,
               &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
               &dim, &p1, &p2, &jj,
               &Dsijj, &Dsij, &Dsjj );

		log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &selected_edge_i, &selected_edge_j,
               &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
               &dim, &p1, &p2, &jj,
               &Dijj, &Dij, &Djj );

		//alpha_ij = ( G[ij] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );
		alpha_ij = ( G[ ij ] ) ? ( logH_ij - logI_p ) - log_ratio_g_prior[ ij ] : ( logI_p - logH_ij ) + log_ratio_g_prior[ ij ];
		// - - - -  End of calculating alpha - - - - - - - - - - - - - - - - - - - - - - - - - - - |
		  		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < alpha_ij )
		{
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

// - - - STEP 3: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );	
		
// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
			for( i = 0; i < pxp ; i++ )
			{
				K_hat_Cpp[ i ]   += K[ i ];
				p_links_Cpp[ i ] += G[ i ];
			}	
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 

	for( i = 0; i < pxp; i++ )
	{	
		p_links[ i ] = p_links_Cpp[ i ];
		K_hat[ i ]   = K_hat_Cpp[ i ];
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian copula Graphical models  
// Based on Double Metropolis-Hastings
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_DMH_rjmcmc_map( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double Ti[], 
                        double K[], int *p, double *threshold, 
                        double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                        int all_graphs[], double all_weights[], double K_hat[], 
                        char *sample_graphs[], double graph_weights[], int *size_sample_g,
                        int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, b1 = *b, count_all_g = 0;
	int index_selected_edge, selected_edge_i, selected_edge_j, counter, size_sample_graph = *size_sample_g;
	int ip, i, j, ij, jj, dim = *p, pxp = dim * dim;
	int qp = dim * ( dim - 1 ) / 2;
	int p1 = dim - 1, p2 = dim - 2, p2x2 = p2 * 2;
	double Dsijj, Dsjj, Dsij, logH_ij, logI_p, Dij, Djj, Dijj, alpha_ij;   

	bool this_one;
	
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<char> char_g( qp );              // char string_g[pp];
	vector<double> K121( 4 ); 																

	vector<double> Kj12( p1 );               // K[j, -j]
	vector<double> sigmaj12( p1 );           // sigma[-j, j]  
	vector<double> sigmaj22( p1 * p1 );      // sigma[-j, -j]
	vector<double> Kj22_inv( p1 * p1 ); 
	vector<double> Kj12xK22_inv( p1 ); 
	vector<double> K12( p2x2 );              // K[e, -e]
	vector<double> sigma11( 4 );             // sigma[e, e]
	vector<double> sigma12( p2x2 );          // sigma[e, -e]
	vector<double> sigma22( p2 * p2 );       // sigma[-e, -e]
	vector<double> sigma11_inv( 4 ); 
	vector<double> sigma21xsigma11_inv( p2x2 ); 
	vector<double> sigma2112( p2 * p2 ); 
	vector<double> K22_inv( p2 * p2 ); 
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
	// - - - - - - - - - - - - - - - - - - - - 
	vector<double> sigma_dmh( pxp );          // for double Metropolis-Hastings
	vector<double> K_dmh( pxp );              // for double Metropolis-Hastings

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		ip = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ ip + j ];
	}

	// For finding the index of rates 
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
	//vector<double> rates( sub_qp );

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			ij = j * dim + i;
			log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
		}

// - - Main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -| 
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
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		// Randomly selecting one edge: NOTE qp = p * ( p - 1 ) / 2 
		index_selected_edge   = static_cast<int>( unif_rand() * sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

		// - - - -  Calculating alpha - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
		ij    = selected_edge_j * dim + selected_edge_i;
		jj    = selected_edge_j * dim + selected_edge_j;
		Dsij  = Ds[ ij ];
		Dsjj  = Ds[ jj ];
		Dsijj = - Dsij * Dsij / Dsjj;
		Dij   = D[ ij ];
		Djj   = D[ jj ];
		Dijj  = - Dij * Dij / Djj;

		// sampling from K and sigma for double Metropolis-Hastings
		rgwish_sigma( G, &size_node[0], Ti, &K_dmh[0], &sigma_dmh[0], &b1, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );		
		
		log_H_ij( &K[0], &sigma[0], &logH_ij, &selected_edge_i, &selected_edge_j,
               &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
               &dim, &p1, &p2, &jj,
               &Dsijj, &Dsij, &Dsjj );

		log_H_ij( &K_dmh[0], &sigma_dmh[0], &logI_p, &selected_edge_i, &selected_edge_j,
               &Kj12[0], &Kj12xK22_inv[0], &K12[0], &K12xK22_inv[0], &K121[0], 
               &sigmaj12[0], &sigmaj22[0], &sigma12[0], &sigma22[0], &sigma11_inv[0], &sigma21xsigma11_inv[0],
               &dim, &p1, &p2, &jj,
               &Dijj, &Dij, &Djj );

		//alpha_ij = ( G[ij] ) ? ( logH_ij - logI_p ) : ( logI_p - logH_ij );
		alpha_ij = ( G[ ij ] ) ? ( logH_ij - logI_p ) - log_ratio_g_prior[ ij ] : ( logI_p - logH_ij ) + log_ratio_g_prior[ ij ];
		// - - - -  End of calculating alpha - - - - - - - - - - - - - - - - - - - - - - - - - - - |

		// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		if( i_mcmc >= burn_in )
		{
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
				all_graphs[ count_all_g ]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
		} 
        // - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	

		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < alpha_ij )
		{
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

// - - - STEP 3: Sampling from G-Wishart for new graph - - - - - - - - - - - - - - - - - - - - - - |
		rgwish_sigma( G, &size_node[0], Ts, K, &sigma[0], b_star, &dim, threshold, &sigma_start[0], &inv_C[0], &beta_star[0], &sigma_i[0], sigma_start_N_i, sigma_N_i, N_i );	
	}  
	PutRNGstate();
// - - End of main loop for birth-death MCMC - - - - - - - - - - - - - - - - - - - - - - - - - - - | 
	
	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
      
} // End of exturn "C"
