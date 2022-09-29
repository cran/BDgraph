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
// Computing the Marginal pseudo-likelihood for HC algorithm and binary data
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_mpl_binary_hc( int *node, int mb_node[], int *size_node, double *log_mpl_node, 
                  int data[], int freq_data[], int *length_freq_data, 
                  double *alpha_ijl, int *n )
{
	double alpha_jl      = 2 * *alpha_ijl;
	double log_alpha_ijl = lgammafn_sign( *alpha_ijl, NULL );
	double log_alpha_jl  = lgammafn_sign( alpha_jl, NULL );

	int i_hash, counter, i, j, l, size_mb_conf, mb_node_x_lf, node_x_lf = *node * *length_freq_data;
	double sum_lgamma_fam;

	*log_mpl_node = 0.0;
	int fam_conf_count_0 = 0, fam_conf_count_1 = 0;
		   
	switch( *size_node )
	{
		case 0:
			for( i = 0; i < *length_freq_data; i++ )
				( data[ node_x_lf + i ] == 0 ) ? fam_conf_count_0 += freq_data[ i ] : fam_conf_count_1 += freq_data[ i ];
				
			sum_lgamma_fam = lgammafn_sign( fam_conf_count_0 + *alpha_ijl, NULL ) + lgammafn_sign( fam_conf_count_1 + *alpha_ijl, NULL );
						   
			*log_mpl_node = sum_lgamma_fam - lgammafn_sign( *n + alpha_jl, NULL ) + log_alpha_jl - 2 * log_alpha_ijl;      
			break;
				
		case 1:
			mb_node_x_lf = mb_node[0] * *length_freq_data; 		

			for( l = 0; l < 2; l++ )  // collects the necessary statistics from the data and calculates the score
			{
				fam_conf_count_0 = 0;
				fam_conf_count_1 = 0;
				for( i = 0; i < *length_freq_data; i++ )
					if( data[ mb_node_x_lf + i ] == l ) 
						( data[ node_x_lf + i ] == 0 ) ? fam_conf_count_0 += freq_data[ i ] : fam_conf_count_1 += freq_data[ i ];
								
				sum_lgamma_fam = lgammafn_sign( fam_conf_count_0 + *alpha_ijl, NULL ) + lgammafn_sign( fam_conf_count_1 + *alpha_ijl, NULL );
		   
				//mb_conf_count = fam_conf_count[0] + fam_conf_count[1];
				*log_mpl_node += sum_lgamma_fam - lgammafn_sign( fam_conf_count_0 + fam_conf_count_1 + alpha_jl, NULL );     
			}		

			// adding remaining terms 
			*log_mpl_node += 2 * log_alpha_jl - 4 * log_alpha_ijl;   
			break;
	
		default:			
			int size_bit = sizeof( unsigned long long int ) * CHAR_BIT / 2;
			int sz = size_bit;
			
			vector<int>vec_fam_conf_count_0( *length_freq_data );	
			vector<int>vec_fam_conf_count_1( *length_freq_data );	

			vector<vector<unsigned long long int> > mb_conf( *length_freq_data );		
			vector<vector<unsigned long long int> > data_mb( *length_freq_data );		

			int size_hash = *size_node / sz + 1;
			vector<unsigned long long int> hash_mb( size_hash, 0 );			

			// for case i = 0
			for( j = 0; j < *size_node; j++ )
			{
				i_hash = j / sz;
				//hash_mb[ i_hash ] |= data[ mb_node[j] * *length_freq_data ] << ( j - i_hash * sz );
				hash_mb[ i_hash ] += (unsigned long long)data[ mb_node[ j ] * *length_freq_data ] << ( j - i_hash * sz );
			}
			mb_conf[0] = hash_mb;
			size_mb_conf = 1;

			if( data[ node_x_lf ] == 0 ) 
			{
				vec_fam_conf_count_0[0] = freq_data[0];
				vec_fam_conf_count_1[0] = 0;
			}else{
				vec_fam_conf_count_1[0] = freq_data[0];
				vec_fam_conf_count_0[0] = 0;
			}

			int sizeof_hash = size_hash * sizeof hash_mb[0];
			for( i = 1; i < *length_freq_data; i++ )
			{
				//vector<unsigned long long int> hash_mb( size_hash, 0 );
				memset( &hash_mb[0], 0, sizeof_hash );
				
				for( j = 0; j < *size_node; j++ )
				{
					i_hash = j / sz;
					//hash_mb[ i_hash ] |= data[ mb_node[j] * *length_freq_data + i ] << ( j - i_hash * sz );
					hash_mb[ i_hash ] += (unsigned long long)data[ mb_node[ j ] * *length_freq_data + i ] << ( j - i_hash * sz );
				}
				
				//data_mb[i] = hash_mb;
				counter = 1;
				for( j = 0; j < size_mb_conf; j++ )
					if( hash_mb == mb_conf[ j ] ) 
					{
						( data[ node_x_lf + i ] == 0 ) ? vec_fam_conf_count_0[ j ] += freq_data[ i ] : vec_fam_conf_count_1[ j ] += freq_data[ i ];
						counter = 0;  
						break;
					}					
				
				if( counter )
				{
					//( data[ node_x_lf + i ] == 0 ) ? vec_fam_conf_count_0[ size_mb_conf ] = freq_data[i] : vec_fam_conf_count_1[ size_mb_conf ] = freq_data[i];
					if( data[ node_x_lf + i ] == 0 ) 
					{
						vec_fam_conf_count_0[ size_mb_conf ] = freq_data[ i ];
						vec_fam_conf_count_1[ size_mb_conf ] = 0;
					}else{
						vec_fam_conf_count_1[ size_mb_conf ] = freq_data[ i ];
						vec_fam_conf_count_0[ size_mb_conf ] = 0;
					}
					
					mb_conf[ size_mb_conf++ ] = hash_mb;
				}
			}
						
			for( l = 0; l < size_mb_conf; l++ )  // collects the necessary statistics from the data and calculates the score
				*log_mpl_node += lgammafn_sign( vec_fam_conf_count_0[l] + *alpha_ijl, NULL ) + lgammafn_sign( vec_fam_conf_count_1[l] + *alpha_ijl, NULL ) - lgammafn_sign( vec_fam_conf_count_0[l] + vec_fam_conf_count_1[l] + alpha_jl, NULL );     

			// adding remaining terms 
			*log_mpl_node += size_mb_conf * ( log_alpha_jl - 2 * log_alpha_ijl );   

			break;
	}
}
    
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Parallel function to compute the Marginal pseudo-likelihood for BINARY data
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
void log_mpl_binary_parallel_hc( int *node, int mb_node[], int *size_node, double *log_mpl_node, 
                  int data[], int freq_data[], int *length_freq_data, 
                  double *alpha_ijl, int *n )
{
	double alpha_jl      = 2 * *alpha_ijl;
	double log_alpha_ijl = lgammafn_sign( *alpha_ijl, NULL );
	double log_alpha_jl  = lgammafn_sign( alpha_jl, NULL );

	int i, l, size_mb_conf, mb_node_x_lf, node_x_lf = *node * *length_freq_data;
	double sum_lgamma_fam;

	*log_mpl_node = 0.0;
	int fam_conf_count_0 = 0, fam_conf_count_1 = 0;
		   
	switch( *size_node )
	{
		case 0:
			for( i = 0; i < *length_freq_data; i++ )
				( data[ node_x_lf + i ] == 0 ) ? fam_conf_count_0 += freq_data[ i ] : fam_conf_count_1 += freq_data[ i ];
				
			sum_lgamma_fam = lgammafn_sign( fam_conf_count_0 + *alpha_ijl, NULL ) + lgammafn_sign( fam_conf_count_1 + *alpha_ijl, NULL );
						   
			*log_mpl_node = sum_lgamma_fam - lgammafn_sign( *n + alpha_jl, NULL ) + log_alpha_jl - 2 * log_alpha_ijl;      
			break;
				
		case 1:
			mb_node_x_lf = mb_node[0] * *length_freq_data; 		

			for( l = 0; l < 2; l++ )  // collects the necessary statistics from the data and calculates the score
			{
				fam_conf_count_0 = 0;
				fam_conf_count_1 = 0;
				for( i = 0; i < *length_freq_data; i++ )
					if( data[ mb_node_x_lf + i ] == l ) 
						( data[ node_x_lf + i ] == 0 ) ? fam_conf_count_0 += freq_data[ i ] : fam_conf_count_1 += freq_data[ i ];
								
				sum_lgamma_fam = lgammafn_sign( fam_conf_count_0 + *alpha_ijl, NULL ) + lgammafn_sign( fam_conf_count_1 + *alpha_ijl, NULL );
		   
				//mb_conf_count = fam_conf_count[0] + fam_conf_count[1];
				*log_mpl_node += sum_lgamma_fam - lgammafn_sign( fam_conf_count_0 + fam_conf_count_1 + alpha_jl, NULL );     
			}		

			// adding remaining terms 
			*log_mpl_node += 2 * log_alpha_jl - 4 * log_alpha_ijl;   
			break;
	
		default:			
			vector<vector<unsigned long long int> > mb_conf( *length_freq_data );		
			vector<vector<unsigned long long int> > data_mb( *length_freq_data );		
			int size_bit = sizeof( unsigned long long int ) * CHAR_BIT / 2;

			int sz = size_bit;
			int size_hash = *size_node / sz + 1;			
	
			#pragma omp parallel
			{
				int j, i_hash;
				vector<unsigned long long int> hash_mb( size_hash );
				int sizeof_hash = size_hash * sizeof hash_mb[0];
				
				#pragma omp for
				for( int i = 0; i < *length_freq_data; i++ )
				{
					memset( &hash_mb[0], 0, sizeof_hash );
					
					for( j = 0; j < *size_node; j++ )
					{
						i_hash = j / sz;
						//hash_mb[ i_hash ] |= data[ mb_node[j] * *length_freq_data + i ] << ( j - i_hash * sz );
						hash_mb[ i_hash ] += (unsigned long long)data[ mb_node[ j ] * *length_freq_data + i ] << ( j - i_hash * sz );
					}
					data_mb[ i ] = hash_mb;
				}
			}
			
			mb_conf = data_mb;
			std::sort( mb_conf.begin(), mb_conf.end() );
			mb_conf.erase( std::unique( mb_conf.begin(), mb_conf.end() ), mb_conf.end() );
			size_mb_conf = mb_conf.size();

			for( l = 0; l < size_mb_conf; l++ )  // collects the necessary statistics from the data and calculates the score
			{
				fam_conf_count_0 = 0;
				fam_conf_count_1 = 0;
				for( i = 0; i < *length_freq_data; i++ )
					if( data_mb[ i ] == mb_conf[ l ] ) 
						( data[ node_x_lf + i ] == 0 ) ? fam_conf_count_0 += freq_data[ i ] : fam_conf_count_1 += freq_data[ i ];
				
				*log_mpl_node += lgammafn_sign( fam_conf_count_0 + *alpha_ijl, NULL ) + lgammafn_sign( fam_conf_count_1 + *alpha_ijl, NULL ) - lgammafn_sign( fam_conf_count_0 + fam_conf_count_1 + alpha_jl, NULL );     
			}

			// adding remaining terms 
			*log_mpl_node += size_mb_conf * ( log_alpha_jl - 2 * log_alpha_ijl );   

			break;
	}
}
    	
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
// Computing the Marginal pseudo-likelihood for discrete data
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
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
		for( j = 0; j < size_mb_conf; j++ ) mb_conf[ j ] = j;     
	}
	
	if( *size_node > 1 ) 
	{
		vector<int>cumprod_mb( *size_node );
		cumprod_mb[0] = max_range_nodes[ mb_node[ 0 ] ];
		//cumprod_mb   = t( t( c( 1, cumprod( max_range_nodes[ mb_node[ 2:length( mb_node ) ] ] ) ) ) )
		for( j = 1; j < *size_node; j++ )
			cumprod_mb[ j ] = cumprod_mb[ j - 1 ] * max_range_nodes[ mb_node[ j ] ];
			
        //data_mb = c( data[ , mb_node ] %*% cumprod_mb ) 
		for( i = 0; i < *length_freq_data; i++ )
			for( j = 0; j < *size_node; j++ )
				data_mb[ i ] += cumprod_mb[ j ] * data[ mb_node[ j ] * *length_freq_data + i ];
			
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
			fam_conf_count[ j ] = 0;
			for( i = 0; i < *length_freq_data; i++ )
				if( data[ node_x_lf + i ] == j ) fam_conf_count[ j ] += freq_data[ i ];
		}
			
		sum_lgamma_fam = 0.0;
		for( j = 0; j < max_range_node_j; j++ ) 
			sum_lgamma_fam += lgammafn( fam_conf_count[ j ] + *alpha_ijl );
					   
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
				if( data_mb[ i ] == mb_conf_l ) mb_conf_count += freq_data[ i ];
			
			//fam_conf_count = node_conf * 0
			//for( j in 1:max_range_node_j ) fam_conf_count[j] = sum( ( data[ , node ] == j ) * ind )
			for( j = 0; j < max_range_node_j; j++ )
			{
				fam_conf_count[ j ] = 0;
				for( i = 0; i < *length_freq_data; i++ )
					if( ( data[ node_x_lf + i ] == j ) && ( data_mb[ i ] == mb_conf_l ) ) fam_conf_count[ j ] += freq_data[ i ];
			}

			sum_lgamma_fam = 0.0;
			for( j = 0; j < max_range_node_j; j++ ) 
				sum_lgamma_fam += lgammafn( fam_conf_count[ j ] + *alpha_ijl );
	   
			*log_mpl_node += sum_lgamma_fam - lgammafn( mb_conf_count + alpha_jl );     
		}		
	}
	
	// adding remaining terms 
    *log_mpl_node += size_mb_conf * lgammafn( alpha_jl ) - size_mb_conf * max_range_node_j * lgammafn( *alpha_ijl );     
}
     
} // End of exturn "C"
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
