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
#     BDMCMC algorithm for graphical models based on marginal pseudo-likelihood
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bdgraph.mpl = function( data, n = NULL, method = "ggm", transfer = TRUE, algorithm = "bdmcmc", 
					iter = 5000, burnin = iter / 2, g.prior = 0.2, g.start = "empty", 
					jump = NULL, alpha = 0.5, save = FALSE, 
					cores = NULL, operator = "or", verbose = TRUE )
{
    if( iter < burnin ) stop( " 'iter' must be higher than 'burnin'" )
    burnin = floor( burnin )

    if( is.numeric( verbose ) )
    {
        if( ( verbose < 1 ) | ( verbose > 100 ) ) 
            stop( "'verbose' (for numeric case) must be between ( 1, 100 )" )
        
        trace_mcmc = floor( verbose )
        verbose = TRUE
    }else{
        trace_mcmc = ifelse( verbose == TRUE, 10, iter + 1000 )
    }
         
	if( inherits( data, "sim" ) ) 
	    data <- data $ data
	
    colnames_data = colnames( data )

	if( !is.matrix( data ) & !is.data.frame( data ) ) 
	    stop( "Data must be a matrix or dataframe" )
	
	if( is.data.frame( data ) ) data <- data.matrix( data )
	
	if( any( is.na( data ) ) ) 
	    stop( "'bdgraph.mpl()' does not deal with missing values. You could use 'bdgraph()' function with option method = 'gcgm'" )	
		
	p <- ncol( data )
	if( p < 3 ) 
	    stop( "Number of variables/nodes ('p') must be more than 2" )
	
	if( is.null( n ) ) 
	    n <- nrow( data )
	
    if( ( is.null( cores ) ) & ( p < 16 ) ) 
        cours = 1
        
    cores = BDgraph::get_cores( cores = cores, verbose = verbose )
	
	if( method == "ggm" ) 
	{
		if( isSymmetric( data ) )
		{
			if ( is.null( n ) ) 
			    stop( "Please specify the number of observations 'n'" )
			
		    cat( "Input is identified as the covariance matrix \n" )
			S <- data
		}else{
 			S <- t( data ) %*% data
		}
	}
   
	if( ( method == "dgm" ) || ( method == "dgm-binary" ) ) 
	{
		if( transfer == TRUE ) 
		    data = transfer( r_data = data )  
	
		p         = ncol( data ) - 1
		freq_data = data[ , p + 1 ]
		data      = data[ , -( p + 1 ) ]
		n         = sum( freq_data )
	
		max_range_nodes = apply( data, 2, max )
		max_range_nodes = max_range_nodes + 1
		length_f_data   = length( freq_data )	
	}
	
	if( method == "dgm-binary" )
		if( ( min( data ) != 0 ) || ( max( data ) != 1 ) ) 
		    stop( "For the case 'method = \"dgm-binary\"', data must be binary, 0 or 1" )
	
	g_prior = BDgraph::get_g_prior( g.prior = g.prior, p = p )
	G       = BDgraph::get_g_start( g.start = g.start, g_prior = g_prior, p = p )
	
	if( save == TRUE )
	{
		qp1           = ( p * ( p - 1 ) / 2 ) + 1
		string_g      = paste( c( rep( 0, qp1 ) ), collapse = '' )
		sample_graphs = c( rep ( string_g, iter - burnin ) )  # vector of numbers like "10100" 
		graph_weights = c( rep ( 1, iter - burnin ) )         # waiting time for every state
		all_graphs    = c( rep ( 0, iter - burnin ) )         # vector of numbers like "10100"
		all_weights   = c( rep ( 1, iter - burnin ) )         # waiting time for every state		
		size_sample_g = 0
	}else{
		p_links = matrix( 0, p, p )
	}

    if( ( verbose == TRUE ) && ( save == TRUE ) && ( p > 50 & iter > 20000 ) )
    {
        cat( "  WARNING: Memory needs to run this function is around: " )
        print( ( iter - burnin ) * utils::object.size( string_g ), units = "auto" ) 
    } 

	last_graph = matrix( 0, p, p )

	if( ( is.null( jump ) ) && ( p > 10 & iter > ( 5000 / p ) ) )
		jump = floor( p / 10 )
	
	if( is.null( jump ) ) 
	    jump = 1
	
	if( ( p < 10 ) && ( jump > 1 ) )      cat( " WARNING: the value of jump should be 1. " )
	if( jump > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: the value of jump should be smaller. " )
	
	if( ( verbose == TRUE ) && ( algorithm != "hc" ) )
		cat( paste( c( iter, " MCMC sampling ... in progress: \n" ), collapse = "" ) ) 

    print = floor( iter / 20 )
	
	# - - - main BDMCMC algorithms implemented in C++ - - - - - - - - - - - - -|
	if( save == TRUE )
	{
		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "ggm_rjmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(S), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}
    
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
		{
			result = .C( "ggm_bdmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(S), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}
    
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
		{
			counter_all_g = 0
			result = .C( "ggm_bdmcmc_mpl_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(S), as.integer(n), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "dgm_rjmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_binary_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
		{
			counter_all_g = 0
			result = .C( "dgm_bdmcmc_mpl_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}
      
		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
		{
			counter_all_g = 0
			result = .C( "dgm_bdmcmc_mpl_binary_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), as.integer(n), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}
      
	}else{
		
		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "ggm_rjmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(S), as.integer(n), as.integer(p), 
						 p_links = as.double(p_links), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
		{
			result = .C( "ggm_bdmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(S), as.integer(n), as.integer(p), 
						 p_links = as.double(p_links), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
		{
			result = .C( "ggm_bdmcmc_mpl_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(S), as.integer(n), as.integer(p), 
						p_links = as.double(p_links), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}				

		if( ( method == "dgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "dgm_rjmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), 
						as.integer(n), as.integer(p), p_links = as.double(p_links), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), 
						as.integer(n), as.integer(p), p_links = as.double(p_links), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_binary_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), 
						as.integer(n), as.integer(p), p_links = as.double(p_links), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p), 
						p_links = as.double(p_links), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}				

		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_binary_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior),  
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), as.integer(n), as.integer(p), 
						p_links = as.double(p_links), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
		}				
	}
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

	if( algorithm != "hc" )
	{
		last_graph             = matrix( result $ G, p, p )
		colnames( last_graph ) = colnames_data[1:p]
   
		if( save == TRUE )
		{
			size_sample_g = result $ size_sample_g
			sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
			graph_weights = result $ graph_weights[ 1 : size_sample_g ]
			all_graphs    = result $ all_graphs + 1
			all_weights   = result $ all_weights	
			
			if( ( algorithm != "rjmcmc" ) & ( jump != 1 ) )
			{ 
				all_weights = all_weights[ 1 : ( result $ counter_all_g ) ]
				all_graphs  = all_graphs[  1 : ( result $ counter_all_g ) ] 
			}

			output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, 
						   all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph,
                           data = data, method = method )
		}else{
			p_links = matrix( result $ p_links, p, p ) 
			
			if( algorithm == "rjmcmc" )	
			    p_links = p_links / ( iter - burnin )
			
			p_links[ lower.tri( p_links ) ] = 0
			colnames( p_links ) = colnames_data[1:p]
			
			output = list( p_links = p_links, last_graph = last_graph,
                           data = data, method = method )
		}
	}else{
		if( method == "dgm" )
			selected_graph = hill_climb_mpl( data = data, freq_data = freq_data, n = n, max_range_nodes = max_range_nodes, alpha = alpha, operator = operator )

		if( method == "dgm-binary" )
			selected_graph = hill_climb_mpl_binary( data = data, freq_data = freq_data, n = n, alpha = alpha, operator = operator )			
		
		colnames( selected_graph ) = colnames_data[ 1:p ]
		output = list( selected_graph = selected_graph,
                       data = data, method = method )
	}
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
	
	class( output ) = "bdgraph"
	return( output )   
}
           
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
