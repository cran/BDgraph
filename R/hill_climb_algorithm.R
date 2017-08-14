## ----------------------------------------------------------------------------|
# A divide-and-conquer type greedy hill climb algorithm
# for undirected graphcial models with dicrete data
# See "Marginal pseudo-likelihood learning of discrete Markov network structures" 
## ----------------------------------------------------------------------------|
# The Hill-Climb algorithm ( function "hill_climb_mpl" ) consists for two part: 
# PART 1: Local Marginal Pseudo-likelihood optimization to discovers the Markov 
#         blanket of each node ( function "local_mb_hc" ).
# PART 2: Neighborhood search algorithm for global Marginal Pseudo-likelihood 
#         optimization  ( function "global_hc" ).
# See "Marginal pseudo-likelihood learning of Markov network structures" by
# Pensar et al. for more details.
## ----------------------------------------------------------------------------|
# INPUT:  * data (n x p) matrix, as discrete data with n observations and p variables. 
#         The outcome space of each variable must be in the form 0, 1, ..., r. 
#         * alpha: The parameter of the prior distribution
# OUTPUT: * selected_g - adjacency matrix for the selected graph 
## ----------------------------------------------------------------------------|
hill_climb_mpl = function( data, freq_data, n, max_range_nodes, alpha = 0.5 )
{
	p = ncol( data )
	
	G = matrix( 0, p, p )	
	for( i in 1:p )
	{
		mes = paste( c( " PART 1: Local search for node ", i ), collapse = "" )
		cat( mes, "\r" )
		flush.console()	
		
		mb_i       = local_mb_hc( node = i, data = data, freq_data = freq_data, max_range_nodes = max_range_nodes, p = p, n = n, alpha = alpha )
		G[mb_i, i] = 1
    }

    G_or = matrix( 0, p, p )
    G_or[ ( G + t( G ) ) > 0 ] = 1    
   
    if( sum( G_or ) != 0 )
    {
       selected_G = global_hc( G_or = G_or, data = data, freq_data = freq_data, max_range_nodes = max_range_nodes, p = p, n = n, alpha = alpha )   
    }else{
		selected_G = G_or
    }
    
    return( selected_G )
}
   
## ----------------------------------------------------------------------------|
# Local Marginal Pseudo-likelihood optimization to discovers the Markov blanket of each node
## ----------------------------------------------------------------------------|
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
   
## ----------------------------------------------------------------------------|
# Neighborhood search algorithm for global Marginal Pseudo-likelihood optimization 
## ----------------------------------------------------------------------------|
global_hc = function( G_or, data, freq_data, max_range_nodes, p, n, alpha = 0.5 )
{
	G_or[ lower.tri( G_or ) ] = 0
	edges_G_or   = which( G_or == 1, arr.ind = T )
	size_G_or    = nrow( edges_G_or )
	length_freq_data = length( freq_data )

	G_hat        = matrix( 0, p, p )
	edges_G_or   = edges_G_or - 1
	
	print( "PART 2, running global search algorithm in C++" )
	
	result = .C( "dgm_global_mpl_hc", G_hat = as.integer(G_hat), as.integer(edges_G_or), as.integer(size_G_or), 
	             as.integer(data), as.integer(freq_data), as.integer(length_freq_data), 
	             as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p), PACKAGE = "BDgraph" )
	
	G_hat = matrix( result $ G_hat, p, p )
	
	return( G_hat )
}
    
## ----------------------------------------------------------------------------|
# Computing the Marginal pseudo-likelihood for discrete data 
## ----------------------------------------------------------------------------|
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
	            as.integer(n), as.integer(p), PACKAGE = "BDgraph" )
	
	log_mpl_node = result $ log_mpl_node
	
	return( log_mpl_node )
}
    









