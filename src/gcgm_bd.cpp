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
#include "rgwish.h"
#include "copula.h"

extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian copula Graphical models  
// for D = I_p 
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_bdmcmc_ma( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double K[], 
            int *p, double *threshold, double Z[], int R[], int not_continuous[], int *n, int *gcgm,
            double K_hat[], double p_links[], int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij;
	int counter = 0, ip, i, j, ij, one = 1, dim = *p, pxp = dim * dim;
	int qp = dim * ( dim - 1 ) / 2;
	double Dsij, weight_C, sum_weights = 0.0, sum_rates; 

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 

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
	vector<double> Dsijj( pxp ); 

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

		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				ij          = j * dim + i;
			    
				Dsij        = Ds[ ij ];
				Dsijj[ ij ] = Dsij * Dsij / Ds[ j * dim + j ]; 
			}
				
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		rates_bdmcmc_parallel( &rates[0], &log_ratio_g_prior[0], G, &index_row[0], &index_col[0], &sub_qp, Ds, &Dsijj[0], &sigma[0], &K[0], b, &dim );

		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			#pragma omp parallel for
			for( i = 0; i < pxp ; i++ )
				if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
			
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
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

	#pragma omp parallel for
	for( i = 0; i < pxp; i++ )
	{	
		p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
		K_hat[ i ]   = K_hat_Cpp[ i ]   / sum_weights;
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian copula Graphical models  
// for D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_bdmcmc_map( int *iter, int *burnin, int G[], double g_prior[], double Ts[], double K[], 
                    int *p, double *threshold, 
                    double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                    int all_graphs[], double all_weights[], double K_hat[], 
                    char *sample_graphs[], double graph_weights[], int *size_sample_g,
                    int *b, int *b_star, double D[], double Ds[], int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, count_all_g = 0;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int counter = 0, ip, i, j, ij, one = 1, dim = *p, pxp = dim * dim;
	bool this_one;
	int qp = dim * ( dim - 1 ) / 2;
	double Dsij, weight_C, sum_weights = 0.0, sum_rates; 

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	vector<char> char_g( qp );              // char string_g[pp];

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

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
	vector<double> Dsijj( pxp ); 

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

		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				ij          = j * dim + i;
				Dsij        = Ds[ ij ];
				Dsijj[ ij ] = Dsij * Dsij / Ds[ j * dim + j ]; 
			}
		
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		rates_bdmcmc_parallel( &rates[0], &log_ratio_g_prior[0], G, &index_row[0], &index_col[0], &sub_qp, Ds, &Dsijj[0], &sigma[0], &K[0], b, &dim );

		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

// - - - Saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			counter = 0;	
			for( j = 1; j < dim; j++ )
				for( i = 0; i < j; i++ )
					char_g[ counter++ ] = G[ j * dim + i ] + '0'; 

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
				all_graphs[ count_all_g ]          = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
		// Updating G (graph) based on selected edge
		selected_edge_ij      = selected_edge_j * dim + selected_edge_i;
		G[ selected_edge_ij ] = 1 - G[ selected_edge_ij];
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

	#pragma omp parallel for
	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;

	#pragma omp parallel for
	for( i = 0; i < pxp; i++ ) 
		K_hat[ i ] /= sum_weights;
}
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Multiple birth-death MCMC for Gaussian copula Graphical models  
// for D = I_p 
// it is for Bayesian model averaging
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_bdmcmc_ma_multi_update( int *iter, int *burnin, int G[], double g_prior[], double Ts[], 
                        double K[], int *p, double *threshold, 
                        double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                        double K_hat[], double p_links[],
                        int *b, int *b_star, double D[], double Ds[], int *multi_update, int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, multi_update_C = *multi_update;
	int selected_edge_i, selected_edge_j, selected_edge_ij;
	int counter = 0, ip, i, j, ij, one = 1, dim = *p, pxp = dim * dim;
	int qp = dim * ( dim - 1 ) / 2;
	double Dsij, weight_C, sum_weights = 0.0, sum_rates; 

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> K_hat_Cpp( pxp, 0.0 ); 
	
	// - - for rgwish_sigma - - - - - - - - - - 
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
	vector<double> Dsijj( pxp ); 

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

		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				ij          = j * dim + i;
				Dsij        = Ds[ ij ];
				Dsijj[ ij ] = Dsij * Dsij / Ds[ j * dim + j ]; 
			}
				
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		rates_bdmcmc_parallel( &rates[0], &log_ratio_g_prior[0], G, &index_row[0], &index_col[0], &sub_qp, Ds, &Dsijj[0], &sigma[0], &K[0], b, &dim );

		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &sub_qp );

// - - - saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			// K_hat_Cpp[i] += K[i] / sum_rates;
			F77_NAME(daxpy)( &pxp, &weight_C, &K[0], &one, &K_hat_Cpp[0], &one );
			
			#pragma omp parallel for
			for( i = 0; i < pxp ; i++ )
				if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
			
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
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

	#pragma omp parallel for
	for( i = 0; i < pxp; i++ )
	{	
		p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
		K_hat[ i ]   = K_hat_Cpp[ i ]   / sum_weights;
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Multiple birth-death MCMC for Gaussian copula Graphical models  
// for D = I_p 
// it is for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void gcgm_bdmcmc_map_multi_update( int *iter, int *burnin, int G[], double g_prior[], double Ts[], 
                    double K[], int *p, double *threshold, 
                    double Z[], int R[], int not_continuous[], int *n, int *gcgm,
                    int all_graphs[], double all_weights[], double K_hat[], 
                    char *sample_graphs[], double graph_weights[], int *size_sample_g, int *counter_all_g,
                    int *b, int *b_star, double D[], double Ds[], int *multi_update, int *print )
{
	int print_c = *print, multi_update_C = *multi_update, iteration = *iter, burn_in = *burnin;
	int count_all_g = *counter_all_g;
	int counter = 0, ip, i, j, ij, one = 1, dim = *p, pxp = dim * dim;
	int selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int qp = dim * ( dim - 1 ) / 2;
	double Dsij, weight_C, sum_weights = 0.0, sum_rates; 
	bool this_one;

	vector<double> sigma( pxp ); 
	vector<double> copyK( pxp ); 
	memcpy( &copyK[0], K, sizeof( double ) * pxp );
	inverse( &copyK[0], &sigma[0], &dim );			

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	vector<char> char_g( qp );              // char string_g[pp];
	
	vector<double> K121( 4 ); 
	// - - for rgwish_sigma - - - - - - - - -
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
	vector<double> Dsijj( pxp ); 

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

		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
			{
				ij          = j * dim + i;
				Dsij        = Ds[ ij ];
				Dsijj[ ij ] = Dsij * Dsij / Ds[ j * dim + j ]; 
			}
				
// - - - STEP 2: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |		

		rates_bdmcmc_parallel( &rates[0], &log_ratio_g_prior[0], G, &index_row[0], &index_col[0], &sub_qp, Ds, &Dsijj[0], &sigma[0], &K[0], b, &dim );

		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &sub_qp );

// - - - Saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |	
		if( i_mcmc >= burn_in )
		{
			counter = 0;	
			for( j = 1; j < dim; j++ )
				for( i = 0; i < j; i++ )
					char_g[ counter++ ] = G[ j * dim + i ] + '0'; 

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

	#pragma omp parallel for
	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy(sample_graphs[ i ], qp, 0);
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	*counter_all_g = count_all_g;

	#pragma omp parallel for
	for( i = 0; i < pxp; i++ ) 
		K_hat[ i ] /= sum_weights;
}
       
} // End of exturn "C"
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
