## ----------------------------------------------------------------------------|
# Main function of BDgraph package: BDMCMC algorithm for graphical models 
## ----------------------------------------------------------------------------|
bdgraph.mpl = function( data, n = NULL, method = "ggm", transfer = TRUE, algorithm = "bdmcmc", 
					iter = 5000, burnin = iter / 2, g.start = "empty", 
					g.space = NULL, g.prior = 0.5, multi.update = NULL, alpha = 0.5, 
					save.all = FALSE, print = 1000, cores = "all", operator = "or" )
{
	if( cores == "all" ) cores = detect_cores()
	.C( "omp_set_num_cores", as.integer( cores ), PACKAGE = "BDgraph" )

	burnin = floor( burnin )
	
	if( class( data ) == "sim" ) data <- data $ data
	colnames_data = colnames( data )

	if( !is.matrix( data ) & !is.data.frame( data ) ) stop( "Data should be a matrix or dataframe" )
	if( is.data.frame( data ) ) data <- data.matrix( data )
	if( iter < burnin ) stop( "Number of iteration must be more than number of burn-in" )
	
	if( ( g.prior <= 0 ) | ( g.prior >= 1 ) ){
		stop( "g.prior must be between 0 and 1" )
	}else{
		g_prior = g.prior
	}

	if( any( is.na( data ) ) ) stop( "This method does not deal with missing values. You could try bdgraph() function with option method = gcgm" )	
		
	p <- ncol( data )
	if( is.null( n ) ) n <- nrow( data )

	if( method == "ggm" ) 
	{
		if( isSymmetric( data ) )
		{
			if ( is.null( n ) ) stop( "Please specify the number of observations 'n'" )
			cat( "Input is identified as the covriance matrix. \n" )
			S <- data
		}else{
 			S <- t( data ) %*% data
		}
	}
   
	if( ( method == "dgm" ) || ( method == "dgm-binary" ) ) 
	{
		if( transfer == TRUE ) data = transfer( r_data = data )  
	
		p         = ncol( data ) - 1
		freq_data = data[ , p + 1 ]
		data      = data[ , -( p + 1 ) ]
		n         = sum( freq_data )
	
		max_range_nodes = apply( data, 2, max )
		max_range_nodes = max_range_nodes + 1
		length_f_data   = length( freq_data )	
	}
	
	if( method == "dgm-binary" )
		if( ( min( data ) != 0 ) || ( max( data ) != 1 ) ) stop( "For the case 'method = dgm-binary', data must be binary (0,1)" )
	
	if( class( g.start ) == "bdgraph"                         ) G = g.start $ last_graph
	if( class( g.start ) == "sim"                             ) G = as.matrix( g.start $ G )
	if( class( g.start ) == "character" && g.start == "empty" ) G = matrix( 0, p, p )
	if( class( g.start ) == "character" && g.start == "full"  )	G = matrix( 1, p, p )
	if( is.matrix( g.start )                                  ) G = g.start
	
	diag(G) = 0
			
	if( save.all == TRUE )
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

	if( ( save.all == TRUE ) && ( p > 50 & iter > 20000 ) )
	{
		cat( "  WARNING: Memory needs to run this function is around " )
		print( ( iter - burnin ) * object.size( string_g ), units = "auto" ) 
	} 
	
	last_graph = matrix( 0, p, p )

	if( ( is.null( multi.update ) ) && ( p > 10 & iter > ( 5000 / p ) ) )
		multi.update = floor( p / 10 )
	
	if( is.null( multi.update ) ) multi.update = 1
	multi_update = multi.update
	
	if( ( p < 10 ) && ( multi_update > 1 ) )      cat( " WARNING: for the cases p < 10, multi.update must be 1 " )
	if( multi_update > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: multi.update must be smaller " )
	
	if( algorithm != "hc" )
	{
		mes <- paste( c( iter, " iteration is started.                    " ), collapse = "" )
		cat( mes, "\r" )
	}
	
	if( is.null( g.space ) )
	{
		g_space = matrix( 1, p, p )
	}else{
		if( !is.matrix( g.space ) ) stop( "g.space must be a matrix (p x p)" )
		g_space = as.matrix( g.space )
	}
	   
## ----------------------------------------------------------------------------|
	if( save.all == TRUE )
	{
		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "ggm_rjmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(S), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(print), PACKAGE = "BDgraph" )
		}
    
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "ggm_bdmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(S), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(print), PACKAGE = "BDgraph" )
		}
    
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "ggm_bdmcmc_mpl_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(S), as.integer(n), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "dgm_rjmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_binary_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), as.integer(n), as.integer(p),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "dgm_bdmcmc_mpl_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}
      
		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "dgm_bdmcmc_mpl_binary_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), as.integer(n), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}
      
	}else{
		
		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "ggm_rjmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(S), as.integer(n), as.integer(p), 
						 p_links = as.double(p_links), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "ggm_bdmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(S), as.integer(n), as.integer(p), 
						 p_links = as.double(p_links), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			result = .C( "ggm_bdmcmc_mpl_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(S), as.integer(n), as.integer(p), 
						p_links = as.double(p_links), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}				

		if( ( method == "dgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "dgm_rjmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), 
						as.integer(n), as.integer(p), p_links = as.double(p_links), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), 
						as.integer(n), as.integer(p), p_links = as.double(p_links), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_binary_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), 
						as.integer(n), as.integer(p), p_links = as.double(p_links), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "dgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), 
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.integer(max_range_nodes), as.double(alpha), as.integer(n), as.integer(p), 
						p_links = as.double(p_links), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}				

		if( ( method == "dgm-binary" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			result = .C( "dgm_bdmcmc_mpl_binary_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior),  
			            as.integer(data), as.integer(freq_data), as.integer(length_f_data), as.double(alpha), as.integer(n), as.integer(p), 
						p_links = as.double(p_links), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}				
	}
## ----------------------------------------------------------------------------|

	if( algorithm != "hc" )
	{
		last_graph             = matrix( result $ G, p, p )
		colnames( last_graph ) = colnames_data[1:p]
   
		if( save.all == TRUE )
		{
			size_sample_g = result $ size_sample_g
			sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
			graph_weights = result $ graph_weights[ 1 : size_sample_g ]
			all_graphs    = result $ all_graphs + 1
			all_weights   = result $ all_weights	
			if( ( algorithm != "rjmcmc" ) & ( multi_update != 1 ) )
			{ 
				all_weights = all_weights[ 1 : ( result $ counter_all_g ) ]
				all_graphs  = all_graphs[  1 : ( result $ counter_all_g ) ] 
			}

			output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, 
						all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph )
		}else{
			p_links = matrix( result $ p_links, p, p ) 
			if( algorithm == "rjmcmc" )	p_links = p_links / ( iter - burnin )
			p_links[ lower.tri( p_links ) ] = 0
			colnames( p_links ) = colnames_data[1:p]
			output = list( p_links = p_links, last_graph = last_graph )
		}
	}else{
		if( method == "dgm" )
			selected_graph = hill_climb_mpl( data = data, freq_data = freq_data, n = n, max_range_nodes = max_range_nodes, alpha = alpha, operator = operator )

		if( method == "dgm-binary" )
			selected_graph = hill_climb_mpl_binary( data = data, freq_data = freq_data, n = n, alpha = alpha, operator = operator )			
		
		colnames( selected_graph ) = colnames_data[1:p]
		output = selected_graph
	}
## ----------------------------------------------------------------------------|
	
	class( output ) = "bdgraph"
	return( output )   
}
           
