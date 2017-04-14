// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2016 Reza Mohammadi
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
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <math.h>        // isinf, sqrt
#include <limits>        // for numeric_limits<double>::max()
#include <R.h>
#include <Rmath.h>
#include <algorithm>     // for transform function and std::sort
#include <functional>    // for transform function

#include "matrix.h"

using namespace std;

extern "C" {
// ----------------------------------------------------------------------------|
// Computing the Marginal pseudo-likelihood for discrete data
// ----------------------------------------------------------------------------|
void log_mpl_dis( int *node, int mb_node[], int *size_node, double *log_mpl_node, 
                  int data[], int freq_data[], int *length_freq_data, 
                  int max_range_nodes[], double *alpha_ijl, int *n )
{
	int i, j, l, size_mb_conf;
	int mb_node_x_lf, mb_conf_l, mb_conf_count;
	int node_x_lf = *node * *length_freq_data;
	double sum_lgamma_fam;

	vector<int>mb_conf( *length_freq_data );
	
    int max_range_node_j = max_range_nodes[ *node ];
	vector<int>fam_conf_count( max_range_node_j );
		   
    double alpha_jl = max_range_node_j * *alpha_ijl;   

	if( *size_node == 0 ) size_mb_conf = 1;
	
	vector<int>data_mb( *length_freq_data, 0 ); 
	
	if( *size_node == 1 ) 
	{
        // data_mb = data[ , mb_node ]          
		mb_node_x_lf = mb_node[0] * *length_freq_data; 
		//for( i = 0; i < *length_freq_data; i++ ) data_mb[i] = data[ mb_node_x_lf + i ]; 
		memcpy( &data_mb[0], &data[0] + mb_node_x_lf, sizeof( int ) * *length_freq_data );    
        
        size_mb_conf = max_range_nodes[ mb_node[0] ];            
        // mb_conf = 1:size_mb_conf;      
		for( j = 0; j < size_mb_conf; j++ ) mb_conf[j] = j + 1;     
	}
	
	if( *size_node > 1 ) 
	{
		vector<int>cumprod_mb( *size_node );
		cumprod_mb[0] = max_range_nodes[ mb_node[ 0 ] ];
		//cumprod_mb   = t( t( c( 1, cumprod( max_range_nodes[ mb_node[ 2:length( mb_node ) ] ] ) ) ) )
		for( j = 1; j < *size_node; j++ )
			cumprod_mb[j] = cumprod_mb[ j - 1 ] * max_range_nodes[ mb_node[ j ] ];
			
        //data_mb = c( data[ , mb_node ] %*% cumprod_mb ) 
		for( i = 0; i < *length_freq_data; i++ )
			for( j = 0; j < *size_node; j++ )
				data_mb[i] += cumprod_mb[j] * data[ mb_node[j] * *length_freq_data + i ];
			
        //mb_conf      = unique( data_mb )            
		//vector<int>mb_conf( *n );
		//for( i = 0; i < *length_freq_data; i++ ) mb_conf[i] = data_mb[i];
		memcpy( &mb_conf[0], &data_mb[0], sizeof( int ) * *length_freq_data );

		std::sort( mb_conf.begin(), mb_conf.end() );
		mb_conf.erase( std::unique( mb_conf.begin(), mb_conf.end() ), mb_conf.end() );	
		
		size_mb_conf = mb_conf.size(); 
	}
	
	if( *size_node == 0 ) 
	{
		//for( i in 1:max_range_node_j ) fam_conf_count[i] = sum( ( data[ , node ] == i ) * ind )
		for( j = 0; j < max_range_node_j; j++ ) 
		{  
			 //fam_conf_count[j] = std::count( &data[0] + node_x_n, &data[0] + node_x_n + *n, j + 1 );
			fam_conf_count[j] = 0;
			for( i = 0; i < *length_freq_data; i++ )
				if( data[node_x_lf + i] == ( j + 1 ) ) fam_conf_count[j] += freq_data[i];
		}
			
		sum_lgamma_fam = 0.0;
		for( j = 0; j < max_range_node_j; j++ ) 
			sum_lgamma_fam += lgammafn( fam_conf_count[j] + *alpha_ijl );
					   
		*log_mpl_node = sum_lgamma_fam - lgammafn( *n + alpha_jl );      
	}

	if( *size_node > 0 ) 
	{
		*log_mpl_node = 0.0;
		for( l = 0; l < size_mb_conf; l++ )  // collects the necessary statistics from the data and calculates the score
		{
			//ind = c( ( data_mb == mb_conf[ l ] ) * 1 )   # finds positions of MB configuration 
			//mb_conf_count = sum( ind )  # n_jl
			mb_conf_l = mb_conf[ l ];
			//mb_conf_count = std::count( data_mb.begin(), data_mb.end(), mb_conf_l );
			mb_conf_count = 0;
			for( i = 0; i < *length_freq_data; i++ )
				if( data_mb[i] == mb_conf_l ) mb_conf_count += freq_data[i];
			
			//fam_conf_count = node_conf * 0
			//for( j in 1:max_range_node_j ) fam_conf_count[j] = sum( ( data[ , node ] == j ) * ind )
			for( j = 0; j < max_range_node_j; j++ )
			{
				fam_conf_count[j] = 0;
				for( i = 0; i < *length_freq_data; i++ )
					if( ( data[node_x_lf + i] == ( j + 1 ) ) && ( data_mb[i] == mb_conf_l ) ) fam_conf_count[j] += freq_data[i];
			}

			sum_lgamma_fam = 0.0;
			for( j = 0; j < max_range_node_j; j++ ) 
				sum_lgamma_fam += lgammafn( fam_conf_count[j] + *alpha_ijl );
	   
			*log_mpl_node += sum_lgamma_fam - lgammafn( mb_conf_count + alpha_jl );     
		}		
	}
	
	// adding remaining terms 
    *log_mpl_node += size_mb_conf * lgammafn( alpha_jl ) - size_mb_conf * max_range_node_j * lgammafn( *alpha_ijl );     
}
    
// ----------------------------------------------------------------------------|
// Computing birth-death rates for dgm_mpl method
// ----------------------------------------------------------------------------|
void rates_gm_mpl_dis( double rates[], double curr_log_mpl[], int G[], int size_node[], int data[], int freq_data[], 
                       int *length_freq_data, int max_range_nodes[], double *alpha_ijl, int *n, int *p )
{
	int i, j, t, nodexdim, counter = 0, count_mb, dim = *p;
	int size_node_i_new, size_node_j_new;
	double log_mpl_i_new, log_mpl_j_new;
	
	double log_rate;

	vector<int>mb_node_i_new( dim );     
	vector<int>mb_node_j_new( dim );     
	
	for( j = 1; j < dim; j++ )
	{
		for( i = 0; i < j; i++ )
		{
			if( G[ j * dim + i ] )
			{ 
				size_node_i_new = size_node[i] - 1; 
				size_node_j_new = size_node[j] - 1; 

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
				size_node_i_new = size_node[i] + 1; 
				size_node_j_new = size_node[j] + 1; 

				nodexdim = i * dim;
				count_mb = 0; 
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] or t == j ) mb_node_i_new[ count_mb++ ] = t;

				nodexdim = j * dim;
				count_mb = 0; 
				for( t = 0; t < dim; t++ ) 
					if( G[ nodexdim + t ] or t == i ) mb_node_j_new[ count_mb++ ] = t;
			}
			
			log_mpl_dis( &i, &mb_node_i_new[0], &size_node_i_new, &log_mpl_i_new, data, freq_data, length_freq_data, max_range_nodes, alpha_ijl, n );		
			log_mpl_dis( &j, &mb_node_j_new[0], &size_node_j_new, &log_mpl_j_new, data, freq_data, length_freq_data, max_range_nodes, alpha_ijl, n );		
																		
			log_rate = log_mpl_i_new + log_mpl_j_new - curr_log_mpl[i] - curr_log_mpl[j];
			rates[ counter++ ] = ( log_rate < 0.0 ) ? exp( log_rate ) : 1.0;			
		}
	}	
}			
   
// ----------------------------------------------------------------------------|
// birth-death MCMC for Graphical models for discrete data with marginal pseudo-likelihood  
// it is for Bayesian model averaging (MA)
// ----------------------------------------------------------------------------|
void dgm_bdmcmc_mpl_ma( int *iter, int *burnin, int G[], int data[], int freq_data[], int *length_f_data, 
               int max_range_nodes[], double *alpha_ijl, int *n, int *p, double p_links[], int *print )
{
	int length_freq_data = *length_f_data;
	int print_c = *print;
	int iteration = *iter, burn_in = *burnin, copy_n = *n;

	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, j, counter;
	int dim = *p, pxp = dim * dim;

	double sum_weights = 0.0, weight_C;
	double sum_rates;
	
	vector<double> p_links_Cpp( pxp, 0.0 ); 
	
	// Counting size of notes
	int nodexdim;
	vector<int>size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ nodexdim + j ];
	}

	// Caclulating the log_likelihood for the current graph G
	int count_mb, t;
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[i] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl_dis( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
	}
	
	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<double> rates( qp );
	vector<int>index_rates_row( qp );
	vector<int>index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

//-- main loop for birth-death MCMC sampling algorithm ------------------------|
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|
				
		rates_gm_mpl_dis( &rates[0], &curr_log_mpl[0], G, &size_node[0], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n, &dim );

//----- Saving result---------------------------- -----------------------------|
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
						
			for( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += weight_C;
			
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	

		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[selected_edge_ij] = 1 - G[selected_edge_ij];
		G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];

		if( G[selected_edge_ij] )
		{ 
			++size_node[selected_edge_i]; 
			++size_node[selected_edge_j]; 
		}
		else
		{ 
			--size_node[selected_edge_i]; 
			--size_node[selected_edge_j]; 
		}

// ----------------------------------------------------------------------------|
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		log_mpl_dis( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		log_mpl_dis( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
// ----------------------------------------------------------------------------|
	} 
	PutRNGstate();
// ----- End of MCMC sampling algorithm ---------------------------------------|

	for( i = 0; i < pxp; i++ )
		p_links[i] = p_links_Cpp[i] / sum_weights;
}
       
// ----------------------------------------------------------------------------|
// birth-death MCMC for Graphical models for discrete data with marginal pseudo-likelihood  
// it is for maximum a posterior probability estimation (MAP)
// ----------------------------------------------------------------------------|
void dgm_bdmcmc_mpl_map( int *iter, int *burnin, int G[], int data[], int freq_data[], int *length_f_data, int max_range_nodes[], double *alpha_ijl, int *n, int *p, 
			 int all_graphs[], double all_weights[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g, int *print )
{
	int length_freq_data = *length_f_data;
	int print_c = *print;
	int iteration = *iter, burn_in = *burnin, copy_n = *n;
	int count_all_g = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	bool this_one;

	int index_selected_edge, selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j, counter;
	int dim = *p;

	double sum_weights = 0.0, weight_C;
	double sum_rates;
	
	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	int nodexdim;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ nodexdim + j ];
	}
	
	// Caclulating the log_likelihood for the current graph G
	int count_mb, t;
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[i] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl_dis( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
	}

	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

//-- main loop for birth-death MCMC sampling algorithm ------------------------|
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc++ )
	{
		if( ( i_mcmc + 1 ) % print_c == 0 ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|
				
		rates_gm_mpl_dis( &rates[0], &curr_log_mpl[0], G, &size_node[0], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n, &dim );
				
//----- Saving result ------------------------ --------------------------------|
		counter = 0;	
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
				char_g[counter++] = G[j * dim + i] + '0'; 

		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[count_all_g] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[count_all_g];
					all_graphs[count_all_g] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[count_all_g];
				all_graphs[count_all_g]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
			
		// Selecting an edge based on birth and death rates
		select_edge( &rates[0], &index_selected_edge, &sum_rates, &qp );
		selected_edge_i = index_rates_row[ index_selected_edge ];
		selected_edge_j = index_rates_col[ index_selected_edge ];

		// Updating G (graph) based on selected edge
		selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
		G[selected_edge_ij] = 1 - G[selected_edge_ij];
		G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];

		if( G[selected_edge_ij] )
		{ 
			++size_node[selected_edge_i]; 
			++size_node[selected_edge_j]; 
		}
		else
		{ 
			--size_node[selected_edge_i]; 
			--size_node[selected_edge_j]; 
		}
// ----------------------------------------------------------------------------|
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		log_mpl_dis( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		log_mpl_dis( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
// ----------------------------------------------------------------------------|
	} 
	PutRNGstate();
// ----- End of MCMC sampling algorithm ---------------------------------------|

	for( i = 0; i < size_sample_graph; i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
}
        
// ----------------------------------------------------------------------------|
// birth-death MCMC for Graphical models for discrete data with marginal pseudo-likelihood  
// it is for Bayesian model averaging (MA)
// ----------------------------------------------------------------------------|
void dgm_bdmcmc_mpl_ma_multi_update( int *iter, int *burnin, int G[], int data[], int freq_data[], int *length_f_data, int max_range_nodes[], double *alpha_ijl, int *n, int *p, 
			 double p_links[], int *multi_update, int *print )
{
	int length_freq_data = *length_f_data;
	int print_c = *print;
	int iteration = *iter, burn_in = *burnin, copy_n = *n;
	int multi_update_C = *multi_update;

	int selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, j, counter;
	int dim = *p, pxp = dim * dim;

	double sum_weights = 0.0, weight_C;
	double sum_rates;

	vector<double> p_links_Cpp( pxp, 0.0 ); 

	// Count size of notes
	int nodexdim;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ nodexdim + j ];
	}

	// Caclulating the log_likelihood for the current graph G
	int count_mb, t;
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[i] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl_dis( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
	}

	// For finding the index of rates 
	int qp = dim * ( dim - 1 ) / 2;
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

//-- main loop for birth-death MCMC sampling algorithm ------------------------|
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % print_c < multi_update_C ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|
				
		rates_gm_mpl_dis( &rates[0], &curr_log_mpl[0], G, &size_node[0], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n, &dim );

//----- Saving result---------------------------- -----------------------------|
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
						
			for( i = 0; i < pxp ; i++ )
				if( G[i] ) p_links_Cpp[i] += weight_C;
			
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	

		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}
// ----------------------------------------------------------------------------|
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		log_mpl_dis( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		log_mpl_dis( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
// ----------------------------------------------------------------------------|
	} 
	PutRNGstate();
// ----- End of MCMC sampling algorithm ---------------------------------------|

	for( i = 0; i < pxp; i++ ) p_links[i] = p_links_Cpp[i] / sum_weights;
}
    
// ----------------------------------------------------------------------------|
// birth-death MCMC for Graphical models for discrete data with marginal pseudo-likelihood  
// it is for maximum a posterior probability estimation (MAP)
// ----------------------------------------------------------------------------|
void dgm_bdmcmc_mpl_map_multi_update( int *iter, int *burnin, int G[], int data[], int freq_data[], int *length_f_data, int max_range_nodes[], double *alpha_ijl, int *n, int *p, 
			 int all_graphs[], double all_weights[], 
			 char *sample_graphs[], double graph_weights[], int *size_sample_g, int *counter_all_g,
			 int *multi_update , int *print )
{
	int length_freq_data = *length_f_data;
	int print_c = *print;
	int multi_update_C = *multi_update;
	int iteration = *iter, burn_in = *burnin, copy_n = *n;
	int count_all_g = *counter_all_g;
	//int count_all_g = 0;
	string string_g;
	vector<string> sample_graphs_C( iteration - burn_in );
	
	bool this_one;

	int selected_edge_i, selected_edge_j, selected_edge_ij, size_sample_graph = *size_sample_g;
	int i, j, counter;
	int dim = *p;

	double sum_weights = 0.0, weight_C;
	double sum_rates;

	int qp = dim * ( dim - 1 ) / 2;
	vector<char> char_g( qp );              // char string_g[pp];

	// Counting size of notes
	int nodexdim;
	vector<int> size_node( dim, 0 );
	for( i = 0; i < dim; i++ )
	{
		nodexdim = i * dim;
		for( j = 0; j < dim; j++ ) size_node[i] += G[ nodexdim + j ];
	}

	// Caclulating the log_likelihood for the current graph G
	int count_mb, t;
	vector<int>mb_node( dim );     
	vector<double>curr_log_mpl( dim );
	for( i = 0; i < dim; i++ ) 
	{ 
		if( size_node[i] > 0 )
		{	
			nodexdim = i * dim;
			count_mb = 0;   
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}
				
		log_mpl_dis( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
	}
	
	// For finding the index of rates 
	vector<double> rates( qp );
	vector<int> index_rates_row( qp );
	vector<int> index_rates_col( qp );
	counter = 0;
	for( j = 1; j < dim; j++ )
		for( i = 0; i < j; i++ )
		{
			index_rates_row[counter] = i;
			index_rates_col[counter] = j;
			counter++;
		}

	int size_index = multi_update_C;
	vector<int> index_selected_edges( multi_update_C );

//-- main loop for birth-death MCMC sampling algorithm ------------------------|
	GetRNGstate();
	for( int i_mcmc = 0; i_mcmc < iteration; i_mcmc += size_index )
	{
		if( ( i_mcmc + 1 ) % print_c < multi_update_C ) Rprintf( " Iteration  %d                 \n", i_mcmc + 1 ); 
		
//----- STEP 1: calculating birth and death rates -----------------------------|
				
		rates_gm_mpl_dis( &rates[0], &curr_log_mpl[0], G, &size_node[0], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n, &dim );
				
//----- Saving result ------------------------ --------------------------------|
		counter = 0;	
		for( j = 1; j < dim; j++ )
			for( i = 0; i < j; i++ )
				char_g[counter++] = G[j * dim + i] + '0'; 
   
		if( i_mcmc >= burn_in )
		{
			weight_C = 1.0 / sum_rates;
			
			string_g = string( char_g.begin(), char_g.end() );	
			all_weights[count_all_g] = weight_C;
			
			this_one = false;
			for( i = 0; i < size_sample_graph; i++ )
				if( sample_graphs_C[i] == string_g )
				{
					graph_weights[i] += all_weights[count_all_g];
					all_graphs[count_all_g] = i;
					this_one = true;
					break;
				} 
			
			if( !this_one || size_sample_graph == 0 )
			{
				sample_graphs_C[size_sample_graph] = string_g;
				graph_weights[size_sample_graph]   = all_weights[count_all_g];
				all_graphs[count_all_g]            = size_sample_graph; 
				size_sample_graph++;				
			}
			
			count_all_g++; 
			sum_weights += weight_C;
		} 
//----- End of saving result --------------------------------------------------|	
			
		// Selecting multiple edges based on birth and death rates
		select_multi_edges( &rates[0], &index_selected_edges[0], &size_index, &sum_rates, &multi_update_C, &qp );

		// Updating graph based on selected edges
		for ( i = 0; i < size_index; i++ )
		{
			selected_edge_i = index_rates_row[ index_selected_edges[i] ];
			selected_edge_j = index_rates_col[ index_selected_edges[i] ];
			
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G[selected_edge_ij] = 1 - G[selected_edge_ij];
			G[selected_edge_i * dim + selected_edge_j] = G[selected_edge_ij];
		
			if( G[selected_edge_ij] )
			{ 
				++size_node[selected_edge_i]; 
				++size_node[selected_edge_j]; 
			}
			else
			{ 
				--size_node[selected_edge_i]; 
				--size_node[selected_edge_j]; 
			}		
		}
// ----------------------------------------------------------------------------|
		if( size_node[ selected_edge_i ] > 0 )
		{	
			nodexdim = selected_edge_i * dim;
			count_mb = 0;  
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}				
		log_mpl_dis( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
		if( size_node[ selected_edge_j ] > 0 )
		{	
			nodexdim = selected_edge_j * dim;
			count_mb = 0;    
			for( t = 0; t < dim; t++ ) 
				if( G[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
		}	
		log_mpl_dis( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
		
// ----------------------------------------------------------------------------|
	} 
	PutRNGstate();
// ----- End of MCMC sampling algorithm ---------------------------------------|

	for( i = 0; i < ( iteration - burn_in ); i++ ) 
	{
		sample_graphs_C[i].copy( sample_graphs[i], qp, 0 );
		sample_graphs[i][qp] = '\0';
	}
	
	*size_sample_g = size_sample_graph;
	*counter_all_g = count_all_g;
}
               
} // End of exturn "C"
