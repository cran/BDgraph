## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2022  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     A divide-and-conquer type greedy hill climb algorithm                    |
#     for undirected graphcial models and count data.                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    The Hill-Climb algorithm (function "hill_climb_mpl") consists for two part:   
#     PART 1: Local Marginal Pseudo-likelihood optimization to discovers the 
#             Markov blanket of each node ( function "local_mb_hc" ).                
#     PART 2: Neighborhood search algorithm for global Marginal Pseudo-likelihood  
#             optimization  ( function "global_hc" ).                               
#     See "Marginal pseudo-likelihood learning of Markov network structures" by     
#     Pensar et al. for more details.                                              
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     INPUT:  * data ( n x p ) matrix, as a count dataset with n observations and p variables. 
#               The outcome space of each variable must be in the form 0, 1, ..., r.   
#             * alpha: The parameter of the prior distribution                        
#     OUTPUT: * selected_G - adjacency matrix for the selected graph                  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

hill_climb_mpl = function( data, freq_data, n, max_range_nodes, alpha = 0.5, operator = "or" )
{
	p = ncol( data )
	
	G = matrix( 0, p, p )	
	for( i in 1:p )
	{
		mes = paste( c( " PART 1: Local search for node ", i ), collapse = "" )
		cat( mes, "\r" )
		utils::flush.console()	
		
		mb_i       = local_mb_hc( node = i, data = data, freq_data = freq_data, max_range_nodes = max_range_nodes, p = p, n = n, alpha = alpha )
		G[mb_i, i] = 1
    }

    G_local = matrix( 0, p, p )
    if( operator == "or" )  G_local[ ( G + t( G ) ) > 0 ] = 1
    if( operator == "and" ) G_local = G * t( G )
  
    if( sum( G_local ) != 0 )
    {
       selected_G = global_hc( G_local = G_local, data = data, freq_data = freq_data, max_range_nodes = max_range_nodes, p = p, n = n, alpha = alpha )   
    }else{
		selected_G = G_local
    }
    
    return( selected_G )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Local Marginal Pseudo-likelihood optimization to discovers the Markov 
#    blanket of each node
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
local_mb_hc = function( node, data, freq_data, max_range_nodes, p, n, alpha = 0.5 )
{
	temp            = seq_len( p )
	mb_potential    = temp[ - node ]	
	l_mb_potential  = p - 1
	mb_hat          = numeric()
	log_prob_mb_hat = log_mpl_disrete( node, mb_hat, data, freq_data, max_range_nodes, p, n, alpha = alpha )	
	cont            = TRUE
	
	while( cont == TRUE )
	{
		cont = FALSE
		log_prob_mb_candidates = numeric( l_mb_potential ) 
		
		for( i in seq_len( l_mb_potential ) )
			log_prob_mb_candidates[i] = log_mpl_disrete( node = node, mb_node = c( mb_hat, mb_potential[i] ), data = data, freq_data = freq_data, max_range_nodes = max_range_nodes, p = p, n = p, alpha = alpha )			
		
		log_prob_mb_candidate_top = max( log_prob_mb_candidates )
		
		if( log_prob_mb_candidate_top > log_prob_mb_hat )
		{
			mb_cand_top_loc = which.max( log_prob_mb_candidates )
			mb_hat          = c( mb_hat, mb_potential[ mb_cand_top_loc ] )
			mb_potential    = mb_potential[ - mb_cand_top_loc ]
			l_mb_potential  = l_mb_potential - 1
			log_prob_mb_hat = log_prob_mb_candidate_top
			cont            = TRUE
		}
		
		length_mb_hat = length( mb_hat )
		if( ( length_mb_hat > 2 ) & ( cont == TRUE ) )
		{
			delete = TRUE
			while( delete == TRUE )
			{
				delete = FALSE
				log_prob_mb_candidates = numeric( length_mb_hat )
				for( i in 1:length_mb_hat )
					log_prob_mb_candidates[i] = log_mpl_disrete( node = node, mb_node = mb_hat[ -i ], data = data, freq_data = freq_data, max_range_nodes = max_range_nodes, p = p, n = n, alpha = alpha )
				
				log_prob_mb_candidate_top = max( log_prob_mb_candidates )
				
				if( log_prob_mb_candidate_top > log_prob_mb_hat )
				{
					mb_cand_top_loc = which.max( log_prob_mb_candidates )
					mb_hat          = mb_hat[ - mb_cand_top_loc ]
					length_mb_hat   = length( mb_hat )
					log_prob_mb_hat = log_prob_mb_candidate_top
					if( length_mb_hat > 2 ) delete = TRUE
				}
			}
		}
	}
	
	return( mb_hat )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Neighborhood search algorithm for global Marginal Pseudo-likelihood optimization 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
global_hc = function( G_local, data, freq_data, max_range_nodes, p, n, alpha = 0.5 )
{
	print( "PART 2, running global search algorithm" )

	ug          = matrix( 0, p, p )
	n_edges     = sum( G_local ) / 2
	temp        = which( G_local == 1, arr.ind = T )
	edges       = temp[temp[, 1] < temp[, 2], ]
	curr_scores = numeric( p )
	
	for( i in 1:p )
		curr_scores[ i ] = log_mpl_disrete( i, which( ug[i, ] == 1 ), data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
	
	edge_change_imp = matrix( 0, n_edges, 2 )
	edge_change     = matrix( TRUE, n_edges, 2 )
	cont = TRUE
	while( cont == TRUE )
	{
		cont = FALSE
		edge_change_ind = which( edge_change[, 1 ] == TRUE )
		
		for( i in 1:length( edge_change_ind ) ) 
		{
			edge = edges[ edge_change_ind[ i ], ]
			node = edge[ 1 ]			
			mb   = which( ug[node, ] == 1 )			
			
			if( ug[ edge[ 1 ], edge[ 2 ] ] == 0 )
			{
				swoe1 = curr_scores[ node ]
				mb    = c( mb, edge[ 2 ] )
				swe1  = log_mpl_disrete( node, mb, data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
				edge_change_imp[ edge_change_ind[ i ], 1 ] = swe1 - swoe1;
				edge_change[ edge_change_ind[ i ], 1 ]     = FALSE
			}else{
				swe1  = curr_scores[ node ]
				mb    = mb[ mb != edge[ 2 ] ]
				swoe1 = log_mpl_disrete( node, mb, data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
				edge_change_imp[ edge_change_ind[ i ], 1 ] = swoe1 - swe1;
				edge_change[ edge_change_ind[ i ], 1 ] = FALSE				
			}			
		}
		
		edge_change_ind = which( edge_change[ , 2 ] == 1 )
		for( i in 1:length( edge_change_ind ) ) 
		{
			edge = edges[ edge_change_ind[ i ], ]
			node = edge[ 2 ]
			mb   = which( ug[ node, ] == 1 )
			
			if( ug[ edge[ 1 ], edge[ 2 ] ] == 0 )
			{
				swoe2 = curr_scores[ node ]
				mb    = c( mb, edge[ 1 ] )
				swe2  = log_mpl_disrete( node, mb, data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
				edge_change_imp[ edge_change_ind[ i ], 2 ] = swe2 - swoe2;
				edge_change[ edge_change_ind[ i ], 2 ] = FALSE
			}else{
				swe2  = curr_scores[ node ]
				mb    = mb[ mb  != edge[ 1 ] ]
				swoe2 = log_mpl_disrete( node, mb, data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
				edge_change_imp[ edge_change_ind[ i ], 2 ] = swoe2 - swe2;
				edge_change[ edge_change_ind[ i ], 2 ]     = FALSE				
			}			
		}
		
		imp         = apply( edge_change_imp, 1, sum )
		max_imp     = max( imp )
		max_imp_loc = which.max( imp )		
		
		if( max_imp > 0 )
		{
			edge = edges[max_imp_loc, ]
			if( ug[edge[1], edge[2]] == 0 )
			{
				ug[edge[1], edge[2]] = 1
				ug[edge[2], edge[1]] = 1
			}else{
				ug[edge[1], edge[2]] = 0
				ug[edge[2], edge[1]] = 0				
			}
			
			curr_scores[edge[1]] = log_mpl_disrete( edge[1], which( ug[edge[1], ] == 1 ), data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
			curr_scores[edge[2]] = log_mpl_disrete( edge[2], which( ug[edge[2], ] == 1 ), data = data, freq_data = freq_data, max_range_nodes, p = p, n = n, alpha = alpha )
			edge_change[edges[, 1] == edge[1]|edges[, 1] == edge[2], 1] = TRUE
			edge_change[edges[, 2] == edge[1]|edges[, 2] == edge[2], 2] = TRUE
			cont = TRUE
		}	
	}	
	
	return( ug )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Computing the Marginal pseudo-likelihood for count data 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
log_mpl_disrete = function( node, mb_node, data, freq_data, max_range_nodes, p, n, alpha = 0.5 )
{
	length_freq_data = length( freq_data )
	size_node        = length( mb_node )
	node             = node - 1
	mb_node          = mb_node - 1
	log_mpl_node     = 0.0

	result = .C( "log_mpl_dis", as.integer(node), as.integer(mb_node), as.integer(size_node), 
	            log_mpl_node = as.double(log_mpl_node), as.integer(data), as.integer(freq_data), 
	            as.integer(length_freq_data), as.integer(max_range_nodes), as.double(alpha), 
	            as.integer(n), PACKAGE = "BDgraph" )
	
	log_mpl_node = result $ log_mpl_node
	
	return( log_mpl_node )
}
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |





