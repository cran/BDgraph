// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//    Copyright (C) 2012 - 2022  Reza Mohammadi                                |
//                                                                             |
//    This file is part of BDgraph package.                                    |
//                                                                             |
//    BDgraph is free software: you can redistribute it and/or modify it under |
//    the terms of the GNU General Public License as published by the Free     |
//    Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
//                                                                             |
//    Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#include "matrix.h"

extern "C" {
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing the Marginal pseudo-likelihood
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_mpl( int *node, int mb_node[], int *size_node, double *log_mpl_node, double S[], 
              double S_mb_node[], int *n, int *p )
{
	int size_node_fa = *size_node + 1, dim = *p, dim1 = dim + 1;
	double det_S_mb_node, det_S_fa_node;

	if( *size_node > 0 )
	{	
		// S_mb_node = S[ mb_node, mb_node ]
		sub_matrix_upper( &S[0], &S_mb_node[0], &mb_node[0], size_node, &dim );
		
		if( *size_node > 1 )
			determinant( &S_mb_node[0], &det_S_mb_node, size_node );
		else
			det_S_mb_node = S[ mb_node[0] * dim1 ];

		// fa_node = c( mb_node, node )
		mb_node[ *size_node ] = *node;   
				
		// S_fa_node = S[fa_node, fa_node]
		sub_matrix_upper( &S[0], &S_mb_node[0], &mb_node[0], &size_node_fa, &dim );
		//det_S_fa_node = det( S_fa_node )
		determinant( &S_mb_node[0], &det_S_fa_node, &size_node_fa );

		//*log_mpl_node = lgammafn( 0.5 * ( *n + *size_node ) ) - lgammafn( 0.5 * size_node_fa ) - ( 2 * *size_node + 1 ) * log( *n ) * 0.5 - ( *n - 1 ) * ( log( det_S_fa_node ) - log( det_S_mb_node ) ) * 0.5;
		*log_mpl_node = lgammafn( 0.5 * ( *n + *size_node ) ) - lgammafn( 0.5 * size_node_fa ) - 
					  *size_node * log( static_cast<double>( *n ) ) - ( *n - 1 ) * ( log( det_S_fa_node ) - log( det_S_mb_node ) ) * 0.5;
	}else{
		det_S_fa_node = S[ *node * dim1 ];
		//*log_mpl_node = lgammafn( 0.5 * *n ) - lgammafn( 0.5 ) - log( *n ) * 0.5 - ( *n - 1 ) * ( log( det_S_fa_node ) ) * 0.5;
		*log_mpl_node = lgammafn( 0.5 * *n ) - lgammafn( 0.5 ) - ( *n - 1 ) * ( log( det_S_fa_node ) ) * 0.5;
	}	
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing birth-death rates for all the possible edges for ggm_mpl method
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void rates_ggm_mpl( double rates[], double log_ratio_g_prior[], double curr_log_mpl[], int G[], 
        int index_row[], int index_col[], int *sub_qp, int size_node[], double S[], int *n, int *p )
{
	int dim = *p;

	#pragma omp parallel
	{
		int i, j, ij, t, nodexdim, count_mb, size_node_i_new, size_node_j_new;
		double log_mpl_i_new, log_mpl_j_new, log_rate_ij;
		
		int *mb_node_i_new = new int[ dim ];           // For dynamic memory used
		int *mb_node_j_new = new int[ dim ];           // For dynamic memory used
		double *S_mb_node  = new double[ dim * dim ];  // For dynamic memory used
		
		#pragma omp for
		for( int counter = 0; counter < *sub_qp; counter++ )
		{
			i  = index_row[ counter ];
			j  = index_col[ counter ];
			ij = j * dim + i;

			if( G[ ij ] )
			{ 
				size_node_i_new = size_node[ i ] - 1; 
				size_node_j_new = size_node[ j ] - 1; 

				if( size_node_i_new > 0 )
				{	
					nodexdim = i * dim;
					count_mb = 0; 
					for( t = 0; t < dim; t++ ) 
						if( G[ nodexdim + t ] and t != j ) mb_node_i_new[ count_mb++ ] = t;
				}	
				
				if( size_node_j_new > 0 )
				{						
					nodexdim = j * dim;
					count_mb = 0; 
					for( t = 0; t < dim; t++ ) 
						if( G[ nodexdim + t ] and t != i ) mb_node_j_new[ count_mb++ ] = t;
				}	
			}else{ 
				size_node_i_new = size_node[ i ] + 1; 
				size_node_j_new = size_node[ j ] + 1; 

				nodexdim = i * dim;
				count_mb = 0; 
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] or t == j ) mb_node_i_new[ count_mb++ ] = t;

				nodexdim = j * dim;
				count_mb = 0; 
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] or t == i ) mb_node_j_new[ count_mb++ ] = t;
			}

			log_mpl( &i, mb_node_i_new, &size_node_i_new, &log_mpl_i_new, S, S_mb_node, n, &dim );		
			log_mpl( &j, mb_node_j_new, &size_node_j_new, &log_mpl_j_new, S, S_mb_node, n, &dim );		
																		
			log_rate_ij = log_mpl_i_new + log_mpl_j_new - curr_log_mpl[ i ] - curr_log_mpl[ j ];
			log_rate_ij = ( G[ ij ] ) ? log_rate_ij - log_ratio_g_prior[ ij ] : log_rate_ij + log_ratio_g_prior[ ij ];
			
			rates[ counter ] = ( log_rate_ij < 0.0 ) ? exp( log_rate_ij ) : 1.0;
		}
		
		delete[] mb_node_i_new;
		delete[] mb_node_j_new;
		delete[] S_mb_node;
	}	
}			
     
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing birth-death rates for all the possible edges for ggm_mpl method
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void local_rates_ggm_mpl( double rates[], double log_ratio_g_prior[], int *selected_edge_i, int *selected_edge_j, 
            double curr_log_mpl[], int G[], int index_row[], int index_col[], int *sub_qp, 
            int size_node[], double S[], int *n, int *p )
{
	int dim = *p;

	#pragma omp parallel
	{
		int i, j, ij, t, nodexdim, count_mb, size_node_i_new, size_node_j_new;
		double log_mpl_i_new, log_mpl_j_new, log_rate_ij;
		
		int *mb_node_i_new = new int[ dim ];           // For dynamic memory used
		int *mb_node_j_new = new int[ dim ];           // For dynamic memory used
		double *S_mb_node  = new double[ dim * dim ];  // For dynamic memory used
		
		#pragma omp for
		for( int counter = 0; counter < *sub_qp; counter++ )
		{
			i  = index_row[ counter ];
			j  = index_col[ counter ];
			ij = j * dim + i;

			if( ( i == *selected_edge_i ) or ( j == *selected_edge_j ) )
			{
				if( G[ ij ] )
				{ 
					size_node_i_new = size_node[ i ] - 1; 
					size_node_j_new = size_node[ j ] - 1; 

					if( size_node_i_new > 0 )
					{	
						nodexdim = i * dim;
						count_mb = 0; 
						for( t = 0; t < dim; t++ ) 
							if( G[ nodexdim + t ] and t != j ) mb_node_i_new[ count_mb++ ] = t;
					}	
					
					if( size_node_j_new > 0 )
					{						
						nodexdim = j * dim;
						count_mb = 0; 
						for( t = 0; t < dim; t++ ) 
							if( G[ nodexdim + t ] and t != i ) mb_node_j_new[ count_mb++ ] = t;
					}	
				}else{ 
					size_node_i_new = size_node[ i ] + 1; 
					size_node_j_new = size_node[ j ] + 1; 

					nodexdim = i * dim;
					count_mb = 0; 
					for( t = 0; t < dim; t++ ) 
						if( G[ nodexdim + t ] or t == j ) mb_node_i_new[ count_mb++ ] = t;

					nodexdim = j * dim;
					count_mb = 0; 
					for( t = 0; t < dim; t++ ) 
						if( G[ nodexdim + t ] or t == i ) mb_node_j_new[ count_mb++ ] = t;
				}

				log_mpl( &i, mb_node_i_new, &size_node_i_new, &log_mpl_i_new, S, S_mb_node, n, &dim );		
				log_mpl( &j, mb_node_j_new, &size_node_j_new, &log_mpl_j_new, S, S_mb_node, n, &dim );		
																			
				log_rate_ij = log_mpl_i_new + log_mpl_j_new - curr_log_mpl[ i ] - curr_log_mpl[ j ];
				log_rate_ij = ( G[ ij ] ) ? log_rate_ij - log_ratio_g_prior[ ij ] : log_rate_ij + log_ratio_g_prior[ ij ];
				
				rates[ counter ] = ( log_rate_ij < 0.0 ) ? exp( log_rate_ij ) : 1.0;
			}
		}
		
		delete[] mb_node_i_new;
		delete[] mb_node_j_new;
		delete[] S_mb_node;
	}	
}			
     
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian Graphical models with marginal pseudo-likelihood  
// for Bayesian model averaging (MA)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_bdmcmc_mpl_ma( int *iter, int *burnin, int G[], double g_prior[], 
                        double S[], int *n, int *p, double p_links[] , int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, copy_n = *n;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij;
	int nodexdim, count_mb, t, i, j, ij, counter, dim = *p, pxp = dim * dim;
	double sum_weights = 0.0, weight_C, sum_rates;
	
	vector<double> p_links_Cpp( pxp, 0.0 ); 
	
	vector<double> copyS( pxp ); 
	memcpy( &copyS[0], S, sizeof( double ) * pxp );

	// Counting size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ nodexdim + j ];
	}
	
	// Caclulating the log_likelihood for the current graph G
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	vector<double> S_mb_node( pxp );     // For dynamic memory used
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[ i ] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	}
	
	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = j * dim + i;
	        log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
	    }

	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<int>index_row( qp );
	vector<int>index_col( qp );
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
	// calculating all the birth and death rates 
	rates_ggm_mpl( &rates[0], &log_ratio_g_prior[0], &curr_log_mpl[0], G, &index_row[0], &index_col[0], &sub_qp, &size_node[0], &copyS[0], &copy_n, &dim );

// - - main loop for birth-death MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - |
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + 1 ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	  						
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &sub_qp );
		selected_edge_i = index_row[ index_selected_edge ];
		selected_edge_j = index_col[ index_selected_edge ];

// - - - Saving result- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
						
			#pragma omp parallel for
			for( i = 0; i < pxp ; i++ )
				if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
			
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	

		// Updating G (graph) based on selected edge
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		log_mpl( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
		
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		log_mpl( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], &copyS[0], &S_mb_node[0], &copy_n, &dim );

		// Calculating local birth and death rates 				
		local_rates_ggm_mpl( &rates[0], &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, &curr_log_mpl[0], G, &index_row[0], &index_col[0], &sub_qp, &size_node[0], &copyS[0], &copy_n, &dim );
		
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	} 
	PutRNGstate();
// - - - End of MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	#pragma omp parallel for
	for( i = 0; i < pxp; i++ )
		p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
}
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// birth-death MCMC for Gaussian Graphical models with marginal pseudo-likelihood  
// for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_bdmcmc_mpl_map( int *iter, int *burnin, int G[], double g_prior[], 
                         double S[], int *n, int *p, int all_graphs[], double all_weights[], 
			             char *sample_graphs[], double graph_weights[], int *size_sample_g , 
			             int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, copy_n = *n, count_all_g = 0;
	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int nodexdim, count_mb, t, i, j, ij, counter, dim = *p, pxp = dim * dim;
	double sum_weights = 0.0, weight_C, sum_rates;
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	vector<double> copyS( pxp ); 
	memcpy( &copyS[0], S, sizeof( double ) * pxp );

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ nodexdim + j ];
	}
	
	// Caclulating the log_likelihood for the current graph G
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	vector<double> S_mb_node( pxp );     // For dynamic memory used
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[ i ] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	}

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = j * dim + i;
	        log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
	    }
	    
	// For finding the index of rates 
	vector<int>index_row( qp );
	vector<int>index_col( qp );
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
	
	// calculating all the birth and death rates 
	rates_ggm_mpl( &rates[0], &log_ratio_g_prior[0], &curr_log_mpl[0], G, &index_row[0], &index_col[0], &sub_qp, &size_node[0], &copyS[0], &copy_n, &dim );

// - - main loop for birth-death MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - |
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + 1 ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	  						
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
			
		// Updating G (graph) based on selected edge
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
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		log_mpl( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
		
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		log_mpl( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	
		// Calculating local birth and death rates 				
		local_rates_ggm_mpl( &rates[0], &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, &curr_log_mpl[0], G, &index_row[0], &index_col[0], &sub_qp, &size_node[0], &copyS[0], &copy_n, &dim );
		
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	} 
	PutRNGstate();
// - - - End of MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	#pragma omp parallel for
	for( i = 0; i < size_sample_graph; i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
        
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Multiple birth-death MCMC for Gaussian Graphical models with marginal pseudo-likelihood  
// for Bayesian model averaging (MA)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_bdmcmc_mpl_ma_multi_update( int *iter, int *burnin, int G[], double g_prior[], 
                        double S[], int *n, int *p, double p_links[], int *multi_update, int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, copy_n = *n;
	int multi_update_C = *multi_update, selected_edge_i, selected_edge_j, selected_edge_ij;
	int nodexdim, count_mb, t, i, j, ij, counter, dim = *p, pxp = dim * dim;
	double sum_weights = 0.0, weight_C, sum_rates;

	vector<double> p_links_Cpp( pxp, 0.0 ); 
	vector<double> copyS( pxp ); 
	memcpy( &copyS[0], S, sizeof( double ) * pxp );

	// Count size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ nodexdim + j ];
	}

	// Caclulating the log_likelihood for the current graph G
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	vector<double> S_mb_node( pxp );     // For dynamic memory used
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[ i ] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	}

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = j * dim + i;
	        log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
	    }
	    
	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<int>index_row( qp );
	vector<int>index_col( qp );
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

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

// - - main loop for birth-death MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - |
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + size_index ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + size_index ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	 				
// - - - STEP 1: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |
				
		rates_ggm_mpl( &rates[0], &log_ratio_g_prior[0], &curr_log_mpl[0], G, &index_row[0], &index_col[0], &sub_qp, &size_node[0], &copyS[0], &copy_n, &dim );
		
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &sub_qp );

// - - - Saving result- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
						
			#pragma omp parallel for
			for( i = 0; i < pxp ; i++ )
				if( G[ i ] ) p_links_Cpp[ i ] += weight_C;
			
			sum_weights += weight_C;
		} 
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	

// - - - Updating graph based on selected edges - - - - - - - - - - - - - - - - - - - - - - - - - -|
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

		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_row[ index_selected_edges[ i ] ];
			selected_edge_j = index_col[ index_selected_edges[ i ] ];

			//curr_log_mpl[ i ] = log_mpl( node = i, mb_node = which( G[ i, ] == 1 ), size_node = sum( G[ i, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
			if( size_node[ selected_edge_i ] > 0 )
			{	
				nodexdim = selected_edge_i * dim;
				count_mb = 0;  
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
			}				
			//log_mpl_dis( &size_node[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_i ], &selected_edge_i, &copy_n, &dim );	
			log_mpl( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
			
			//curr_log_mpl[ j ] = log_mpl( node = j, mb_node = which( G[ j, ] == 1 ), size_node = sum( G[ j, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
			if( size_node[ selected_edge_j ] > 0 )
			{	
				nodexdim = selected_edge_j * dim;
				count_mb = 0;    
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
			}	
			//log_mpl_dis( &size_node[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_j ], &selected_edge_j, &copy_n, &dim );	
			log_mpl( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
		}
		
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	} 
	PutRNGstate();
// - - - End of MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	#pragma omp parallel for
	for( i = 0; i < pxp; i++ ) 
		p_links[ i ] = p_links_Cpp[ i ] / sum_weights;
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Multiple birth-death MCMC for Gaussian Graphical models with marginal pseudo-likelihood  
// for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_bdmcmc_mpl_map_multi_update( int *iter, int *burnin, int G[], double g_prior[], 
                double S[], int *n, int *p, int all_graphs[], double all_weights[], 
                char *sample_graphs[], double graph_weights[], int *size_sample_g, int *counter_all_g,
                int *multi_update , int *print )
{
	int print_c = *print, multi_update_C = *multi_update;
	int iteration = *iter, burn_in = *burnin, copy_n = *n, count_all_g = *counter_all_g;
	int selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int nodexdim, count_mb, t, i, j, ij, counter, dim = *p, pxp = dim * dim;
	double sum_weights = 0.0, weight_C, sum_rates;
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	vector<double> copyS( pxp ); 
	memcpy( &copyS[0], S, sizeof( double ) * pxp );

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ nodexdim + j ];
	}
	
	// Caclulating the log_likelihood for the current graph G
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	vector<double> S_mb_node( pxp );     // For dynamic memory used
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[ i ] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	}

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = j * dim + i;
	        log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
	    }
	    
	// For finding the index of rates 
	vector<int>index_row( qp );
	vector<int>index_col( qp );
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

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

// - - main loop for birth-death MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - |
	GetRNGstate();
	int print_conter = 0;
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + size_index ) % ( ( print_c * iteration ) / 100 ) == 0 ){
		    ++print_conter;
		    ( ( i_mcmc + size_index ) != iteration ) ? Rprintf( "%i%%->", print_c * print_conter ) : Rprintf( " done" );
		}
	 				
// - - - STEP 1: calculating birth and death rates - - - - - - - - - - - - - - - - - - - - - - - - |
				
		rates_ggm_mpl( &rates[0], &log_ratio_g_prior[0], &curr_log_mpl[0], G, &index_row[0], &index_col[0], &sub_qp, &size_node[0], &copyS[0], &copy_n, &dim );
		
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

// - - - Updating graph based on selected edges - - - - - - - - - - - - - - - - - - - - - - - - - -|
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

		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_row[ index_selected_edges[ i ] ];
			selected_edge_j = index_col[ index_selected_edges[ i ] ];

			//curr_log_mpl[ i ] = log_mpl( node = i, mb_node = which( G[ i, ] == 1 ), size_node = sum( G[ i, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
			if( size_node[ selected_edge_i ] > 0 )
			{	
				nodexdim = selected_edge_i * dim;
				count_mb = 0;  
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
			}				
			//log_mpl_dis( &size_node[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_i ], &selected_edge_i, &copy_n, &dim );	
			log_mpl( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
			
			//curr_log_mpl[ j ] = log_mpl( node = j, mb_node = which( G[ j, ] == 1 ), size_node = sum( G[ j, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
			if( size_node[ selected_edge_j ] > 0 )
			{	
				nodexdim = selected_edge_j * dim;
				count_mb = 0;    
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
			}	
			//log_mpl_dis( &size_node[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_j ], &selected_edge_j, &copy_n, &dim );	
			log_mpl( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
		}
		
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	} 
	PutRNGstate();
// - - - End of MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	#pragma omp parallel for
	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	*counter_all_g = count_all_g;		
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing alpha (probability of acceptness) in RJ-MCMC algorithm for ggm_mpl method
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_alpha_rjmcmc_ggm_mpl( double *log_alpha_ij, double log_ratio_g_prior[], int *i, int *j, 
                        double curr_log_mpl[], int G[], int size_node[], double S[], int *n, int *p )
{
	int t, nodexdim, dim = *p, count_mb, size_node_i_new, size_node_j_new;
	double log_mpl_i_new, log_mpl_j_new;

	vector<int> mb_node_i_new( dim );       // For dynamic memory used
	vector<int> mb_node_j_new( dim );       // For dynamic memory used
	vector<double> S_mb_node( dim * dim );  // For dynamic memory used
	
	int ij = *j * dim + *i;
	if( G[ ij ] )
	{ 
		size_node_i_new = size_node[ *i ] - 1; 
		size_node_j_new = size_node[ *j ] - 1; 

		if( size_node_i_new > 0 )
		{	
			nodexdim = *i * dim;
			count_mb = 0; 
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] and t != *j ) mb_node_i_new[ count_mb++ ] = t;
		}	
		
		if( size_node_j_new > 0 )
		{						
			nodexdim = *j * dim;
			count_mb = 0; 
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] and t != *i ) mb_node_j_new[ count_mb++ ] = t;
		}	
	}else{ 
		size_node_i_new = size_node[ *i ] + 1; 
		size_node_j_new = size_node[ *j ] + 1; 

		nodexdim = *i * dim;
		count_mb = 0; 
		for( t = 0; t < dim; t++ ) 
			if( G[ nodexdim + t ] or t == *j ) mb_node_i_new[ count_mb++ ] = t;

		nodexdim = *j * dim;
		count_mb = 0; 
		for( t = 0; t < dim; t++ ) 
			if( G[ nodexdim + t ] or t == *i ) mb_node_j_new[ count_mb++ ] = t;
	}

	log_mpl( i, &mb_node_i_new[0], &size_node_i_new, &log_mpl_i_new, S, &S_mb_node[0], n, &dim );		
	log_mpl( j, &mb_node_j_new[0], &size_node_j_new, &log_mpl_j_new, S, &S_mb_node[0], n, &dim );		
																
	*log_alpha_ij = log_mpl_i_new + log_mpl_j_new - curr_log_mpl[ *i ] - curr_log_mpl[ *j ];
	*log_alpha_ij = ( G[ ij ] ) ? *log_alpha_ij - log_ratio_g_prior[ ij ] : *log_alpha_ij + log_ratio_g_prior[ ij ];
}			
     
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian Graphical models with marginal pseudo-likelihood  
// for Bayesian model averaging (MA)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_rjmcmc_mpl_ma( int *iter, int *burnin, int G[], double g_prior[], 
                        double S[], int *n, int *p, double p_links[] , int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, copy_n = *n;
	int selected_edge, selected_edge_i, selected_edge_j;
	int nodexdim, count_mb, t, i, j, ij, counter, dim = *p, pxp = dim * dim;
	double log_alpha_ij;
	
	vector<double> p_links_Cpp( pxp, 0.0 ); 
	
	vector<double> copyS( pxp ); 
	memcpy( &copyS[0], S, sizeof( double ) * pxp );

	// Counting size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ nodexdim + j ];
	}
	
	// Caclulating the log_likelihood for the current graph G
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	vector<double> S_mb_node( pxp );     // For dynamic memory used
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[ i ] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	}

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = j * dim + i;
	        log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
	    }
	    
	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<int>index_row( qp );
	vector<int>index_col( qp );
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

// - - main loop for birth-death MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - |
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

		log_alpha_rjmcmc_ggm_mpl( &log_alpha_ij, &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, &curr_log_mpl[0], G, &size_node[0], &copyS[0], &copy_n, &dim );
		
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

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
		//curr_log_mpl[ i ] = log_mpl( node = i, mb_node = which( G[ i, ] == 1 ), size_node = sum( G[ i, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		//log_mpl_dis( &size_node[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_i ], &selected_edge_i, &copy_n, &dim );	
		log_mpl( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
		
		//curr_log_mpl[ j ] = log_mpl( node = j, mb_node = which( G[ j, ] == 1 ), size_node = sum( G[ j, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		//log_mpl_dis( &size_node[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_j ], &selected_edge_j, &copy_n, &dim );	
		log_mpl( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], &copyS[0], &S_mb_node[0], &copy_n, &dim );

// - - - Saving result- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
		if( i_mcmc >= burn_in )
			for( i = 0; i < pxp ; i++ )
				p_links_Cpp[ i ] += G[ i ];
// - - - End of saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|	
		
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
	} 
	PutRNGstate();
// - - - End of MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	memcpy( &p_links[0], &p_links_Cpp[0], sizeof( double ) * pxp );    
}
       
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Reversible Jump MCMC for Gaussian Graphical models with marginal pseudo-likelihood  
// for maximum a posterior probability estimation (MAP)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void ggm_rjmcmc_mpl_map( int *iter, int *burnin, int G[], double g_prior[], double S[], 
                int *n, int *p, int all_graphs[], double all_weights[], 
                char *sample_graphs[], double graph_weights[], int *size_sample_g , int *print )
{
	int print_c = *print, iteration = *iter, burn_in = *burnin, copy_n = *n, count_all_g = 0;
	int selected_edge, selected_edge_i, selected_edge_j, size_sample_graph = *size_sample_g;
	int nodexdim, count_mb, t, i, j, ij, counter, dim = *p, pxp = dim * dim;
	double log_alpha_ij;
	bool this_one;

	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	vector<double> copyS( pxp ); 
	memcpy( &copyS[0], S, sizeof( double ) * pxp );

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[ i ] += G[ nodexdim + j ];
	}
	
	// Caclulating the log_likelihood for the current graph G
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	vector<double> S_mb_node( pxp );     // For dynamic memory used
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[ i ] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], &copyS[0], &S_mb_node[0], &copy_n, &dim );
	}

	vector<double> log_ratio_g_prior( pxp );	
	for( j = 1; j < dim; j++ )
	    for( i = 0; i < j; i++ )
	    {
	        ij = j * dim + i;
	        log_ratio_g_prior[ ij ] = log( static_cast<double>( g_prior[ ij ] / ( 1 - g_prior[ ij ] ) ) );
	    }
	    
	// For finding the index of rates 
	vector<int>index_row( qp );
	vector<int>index_col( qp );
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

// - - main loop for RJ-MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - -|
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

		log_alpha_rjmcmc_ggm_mpl( &log_alpha_ij, &log_ratio_g_prior[0], &selected_edge_i, &selected_edge_j, &curr_log_mpl[0], G, &size_node[0], &copyS[0], &copy_n, &dim );
		
// - - - End of calculating log_alpha_ij - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |		
		  		
		// Selecting an edge and updating G (graph)
		if( log( static_cast<double>( unif_rand() ) ) < log_alpha_ij )
		{
			ij      = selected_edge_j * dim + selected_edge_i;
			G[ ij ] = 1 - G[ ij ];
			G[selected_edge_i * dim + selected_edge_j] = G[ ij ];

			if( G[ ij ] )
			{ 
				++size_node[ selected_edge_i ]; 
				++size_node[ selected_edge_j ]; 
			}else{
				--size_node[ selected_edge_i ]; 
				--size_node[ selected_edge_j ]; 
			}
		}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
		//curr_log_mpl[ i ] = log_mpl( node = i, mb_node = which( G[ i, ] == 1 ), size_node = sum( G[ i, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		//log_mpl_dis( &size_node[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_i ], &selected_edge_i, &copy_n, &dim );	
		log_mpl( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], &copyS[0], &S_mb_node[0], &copy_n, &dim );
		
		//curr_log_mpl[ j ] = log_mpl( node = j, mb_node = which( G[ j, ] == 1 ), size_node = sum( G[ j, ] ), S = S, n = n, p = p, alpha_ijl = alpha_ijl )
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		//log_mpl_dis( &size_node[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, &mb_node[0], alpha_ijl, &curr_log_mpl[ selected_edge_j ], &selected_edge_j, &copy_n, &dim );	
		log_mpl( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], &copyS[0], &S_mb_node[0], &copy_n, &dim );

// - - - Saving result - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
		if( i_mcmc >= burn_in )
		{
			counter = 0;	
			for( j = 1; j < dim; j++ )
				for( i = 0; i < j; i++ )
					char_g[ counter++ ] = G[ j * dim + i ] + '0'; 

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
// - - - End of MCMC sampling algorithm - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	#pragma omp parallel for
	for( i = 0; i < size_sample_graph; i++ ) 
	{
		sample_graphs_C[ i ].copy( sample_graphs[ i ], qp, 0 );
		sample_graphs[ i ][ qp ] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
                 
} // End of exturn "C"
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
