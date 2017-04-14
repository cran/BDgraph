// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2017 Reza Mohammadi
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
#include <algorithm>     // for transform function and std::sort, std::max_element
#include <functional>    // for transform function

#include "matrix.h"

using namespace std;

extern "C" {
// ----------------------------------------------------------------------------|
// Computing the Marginal pseudo-likelihood for discrete data
// ----------------------------------------------------------------------------|
void log_mpl_hc_dis( int *node, int mb_node[], int *size_node, double *log_mpl_node, 
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
// Computing Bayes-factor for Hill-Climb algorithm for discrete data
// ----------------------------------------------------------------------------|
void bayes_factors_mpl_dis( double bayes_factors[], double curr_log_mpl[], int G_hat[], int edges_G_or[], 
                            int *size_G_or, int size_node[], int data[], int freq_data[], 
                            int *length_freq_data, int max_range_nodes[], double *alpha_ijl, int *n, int *p )
{
	int i, j, t, nodexdim, count_mb, dim = *p;
	int size_node_i_new, size_node_j_new;
	double log_mpl_i_new, log_mpl_j_new;
	
	vector<int>mb_node_i_new( dim );     
	vector<int>mb_node_j_new( dim );     
	
	for( int e = 0; e < *size_G_or; e++ )
	{
		i = edges_G_or[ e ];
		j = edges_G_or[ *size_G_or + e ];

		if( G_hat[ j * dim + i ] )
		{ 
			size_node_i_new = size_node[i] - 1; 
			size_node_j_new = size_node[j] - 1; 

			if( size_node_i_new > 0 )
			{	
				nodexdim = i * dim;
				count_mb = 0; 
				for( t = 0; t < dim; t++ ) 
					if( G_hat[ nodexdim + t ] and t != j ) mb_node_i_new[ count_mb++ ] = t;
			}	
			
			if( size_node_j_new > 0 )
			{						
				nodexdim = j * dim;
				count_mb = 0; 
				for( t = 0; t < dim; t++ ) 
					if( G_hat[ nodexdim + t ] and t != i ) mb_node_j_new[ count_mb++ ] = t;
			}	
		}else{ 
			size_node_i_new = size_node[i] + 1; 
			size_node_j_new = size_node[j] + 1; 

			nodexdim = i * dim;
			count_mb = 0; 
			for( t = 0; t < dim; t++ ) 
				if( G_hat[ nodexdim + t ] or t == j ) mb_node_i_new[ count_mb++ ] = t;

			nodexdim = j * dim;
			count_mb = 0; 
			for( t = 0; t < dim; t++ ) 
				if( G_hat[ nodexdim + t ] or t == i ) mb_node_j_new[ count_mb++ ] = t;
		}
		
		log_mpl_hc_dis( &i, &mb_node_i_new[0], &size_node_i_new, &log_mpl_i_new, data, freq_data, length_freq_data, max_range_nodes, alpha_ijl, n );		
		log_mpl_hc_dis( &j, &mb_node_j_new[0], &size_node_j_new, &log_mpl_j_new, data, freq_data, length_freq_data, max_range_nodes, alpha_ijl, n );		
																	
		bayes_factors[ e ] = log_mpl_i_new + log_mpl_j_new - curr_log_mpl[i] - curr_log_mpl[j];
	}	
}			
   
// ----------------------------------------------------------------------------|
// Neighborhood search algorithm for global Marginal Pseudo-likelihood optimization 
// ----------------------------------------------------------------------------|
void dgm_global_mpl_hc( int G_hat[], int edges_G_or[], int *size_G_or, int data[], int freq_data[], int *length_f_data, 
                        int max_range_nodes[], double *alpha_ijl, int *n, int *p )
{
	int length_freq_data = *length_f_data, copy_n = *n;

	int max_bayes_factors_loc, selected_edge_i, selected_edge_j, selected_edge_ij;
	int i, dim = *p;

	// Caclulating the log_likelihood for the current graph G
	int nodexdim, count_mb, t;
	vector<int>size_node( dim, 0 );
	vector<int>mb_node( dim, 0 );     
	vector<double>curr_log_mpl( dim );
	for( i = 0; i < dim; i++ ) 				
		log_mpl_hc_dis( &i, &mb_node[0], &size_node[i], &curr_log_mpl[i], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
	
	vector<double>bayes_factors( *size_G_or );

//-- main loop for neighborhood sampling algorithm ----------------------------|
	int cont = 1;
	while( cont == 1 )
	{
		cont = 0;

		bayes_factors_mpl_dis( &bayes_factors[0], &curr_log_mpl[0], G_hat, edges_G_or, size_G_or, &size_node[0], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n, &dim );
		
		double max_bayes_factors = *std::max_element( bayes_factors.begin(), bayes_factors.end() );
		//double max_bayes_factors = *std::max_element( &bayes_factors[0], &bayes_factors[0] + *size_G_or );
		
		if( max_bayes_factors > 0 )
		{
			//max_bayes_factors_loc = which.max( bayes_factors )
			max_bayes_factors_loc = std::distance( bayes_factors.begin(), std::max_element( bayes_factors.begin(), bayes_factors.end() ) );
			
			selected_edge_i = edges_G_or[ max_bayes_factors_loc ];
			selected_edge_j = edges_G_or[ *size_G_or + max_bayes_factors_loc ];

			// Updating G (graph) based on selected edge
			selected_edge_ij    = selected_edge_j * dim + selected_edge_i;
			G_hat[selected_edge_ij] = 1 - G_hat[selected_edge_ij];
			G_hat[selected_edge_i * dim + selected_edge_j] = G_hat[selected_edge_ij];
						
			if( G_hat[ selected_edge_ij ] )
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
					if( G_hat[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
			}				
			log_mpl_hc_dis( &selected_edge_i, &mb_node[0], &size_node[ selected_edge_i ], &curr_log_mpl[ selected_edge_i ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );
			
			if( size_node[ selected_edge_j ] > 0 )
			{	
				nodexdim = selected_edge_j * dim;
				count_mb = 0;    
				for( t = 0; t < dim; t++ ) 
					if( G_hat[ nodexdim + t ] ) mb_node[ count_mb++ ] = t;
			}	
			log_mpl_hc_dis( &selected_edge_j, &mb_node[0], &size_node[ selected_edge_j ], &curr_log_mpl[ selected_edge_j ], data, freq_data, &length_freq_data, max_range_nodes, alpha_ijl, &copy_n );	
// ----------------------------------------------------------------------------|
			cont = 1;
		}	
// ----------------------------------------------------------------------------|
	} 
// ----- End of sampling algorithm --------------------------------------------|
}
       
} // End of exturn "C"
