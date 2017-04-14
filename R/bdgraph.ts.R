## Main function: BDMCMC algorithm for time series
################################################################################
# Nlength is the length of the time series
# data is the aggregate periodogram Pk, which is arranged as a large p x (Nlength*p) matrix [P1, P2, ... ,PNlength]

bdgraph.ts = function( data, Nlength = NULL, n, iter = 1000, burnin = iter / 2, 
					   g.start = "empty", prior.df = rep( 3, Nlength ), save.all = FALSE )
{
	burnin    = floor( burnin )

	if( !is.matrix(data) & !is.data.frame(data) ) stop( "Data should be a matrix or dataframe" )
	if( is.data.frame(data) ) data <- data.matrix(data)
	if( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )

	if( any( is.na( data ) ) ) stop( "This method does not deal with missing value" )
		
	dimd <- dim(data)
	p    <- dimd[1]
	if( is.null(Nlength) ) Nlength <- dimd[2] / p
	
	P = data 
	W = matrix( 0, p, p * Nlength )  # I(p*p)
	
	for( t in 1 : Nlength )
	  diag( W[, ( t * p - p + 1 ):( t * p )] ) = 1
	
	K      = W
	sigma  = K
	b      = prior.df
	b_star = b + n
	Ws     = W + P
	
	Ls = matrix( 0, 2 * p, 2 * p * Nlength ) # Ls = t(chol([Re(inv_Pk), -Im[inv_Pk]; Re(inv_Pk), Im[inv_Pk]]))
	
	for( k in 1:Nlength )
	{
	  inv_Ws = solve( Ws[, ( k * p - p + 1 ):( k * p )] )
	  row1   = cbind( Re( inv_Ws ), -Im( inv_Ws ) )
	  row2   = cbind( Im( inv_Ws ), Re( inv_Ws ) )
	  sig    = rbind( row1, row2 ) / 2
	  Ls[, ( k * 2 * p - 2 * p + 1 ):( k * 2 * p )] = t( chol( sig ) )
	}

	if( class( g.start ) == "bdgraph" ) 
	{
		G <- g.start $ last_graph
		K <- g.start $ last_K
		if( dim(K)[2] != p * Nlength ) stop( "K should be a p x (p x Nlength) matrix")
	} 

	if( class( g.start ) == "sim" ) 
	{
		G <- as.matrix( g.start $ G )
		K <- as.matrix( g.start $ K )
		if( dim(K)[2] != p * Nlength ) stop( "K should be a p x (p x Nlength) matrix")
	} 
	
	if( class( g.start ) == "character" && g.start == "empty" )
	{
		G = matrix( 0, p, p )
		
		for( t in 1:Nlength )
		{
		  result = .C( "rgcwish_c", as.integer(G), as.double(Ls), K = as.complex(K), as.integer(b_star), as.integer(p), PACKAGE = "BDgraph" )
		  K[, ( t * p - p + 1 ):( t * p )] = matrix ( result $ K, p, p )
		  sigma[, ( t * p - p + 1 ):( t * p )] = solve( K[, ( t * p - p + 1 ):( t * p )] )
		}
	}
	
	if( class( g.start ) == "character" && g.start == "full" )
	{
		G       = matrix( 1, p, p )
		diag(G) = 0
		
		for( t in 1:Nlength )
		{
		  result = .C( "rcwish_c", as.double(Ls), K = as.complex(K), as.integer(b_star), as.integer(p), PACKAGE = "BDgraph" )
		  K[, ( t * p - p + 1 ):( t * p )] = matrix ( result $ K, p, p ) 
		  sigma[, ( t * p - p + 1 ):( t * p )] = solve( K[, ( t * p - p + 1 ):( t * p ) ] )
		}	
	}

	if( is.matrix( g.start ) )
	{
		G       = g.start
		diag(G) = 0
		
		for( t in 1:Nlength )
		{
		  result = .C( "rgcwish_c", as.integer(G), as.double(Ls), K = as.complex(K), as.integer(b_star), as.integer(p), PACKAGE = "BDgraph" )
		  K[, ( t * p - p + 1 ):( t * p )] = matrix ( result $ K, p, p )
		  sigma[, ( t * p - p + 1 ):( t * p )] = solve( K[, ( t * p - p + 1 ):( t * p )] )
		}	
	}

	if( save.all == TRUE )
	{
		qp1           = ( p * ( p - 1 ) / 2 ) + 1
		string_g      = paste( rep( 0, qp1 ), collapse = '' )
		sample_graphs = rep ( string_g, iter - burnin )  # vector of numbers like "10100" 
		graph_weights = rep ( 0, iter - burnin )         # waiting time for every state
		all_graphs    = rep ( 0, iter - burnin )         # vector of numbers like "10100"
		all_weights   = rep ( 1, iter - burnin )         # waiting time for every state		
		size_sample_g = 0
		exit = 0
	}
	else
		p_links = 0 * G

	if( ( save.all == TRUE ) && ( p > 50 & iter > 20000 ) )
	{
		cat( "  WARNING: Memory needs to run this function is around " )
		print( ( iter - burnin ) * object.size( string_g ), units = "auto" ) 
	} 

	r_K        = Re(K)
	i_K        = Im(K)
	r_K_hat    = 0 * K
	i_K_hat    = 0 * K
	r_sigma    = Re(sigma)
	i_sigma    = Im(sigma)
	last_graph = G

	mes <- paste( c( iter, " iteration is started.                    " ), collapse = "" )
	cat( mes, "\r" )

	if ( save.all == TRUE )
	{
		result = .C( "bdmcmc_map_for_multi_dim", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ls), r_K = as.double(r_K), 
					i_K = as.double(i_K), as.integer(p), as.integer(Nlength), r_sigma = as.double(r_sigma), i_sigma = as.double(i_sigma), 
					all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), r_K_hat = as.double(r_K_hat), i_K_hat = as.double(i_K_hat), 
					sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), 
					exit = as.integer(exit), as.integer(b), as.integer(b_star), r_Ds = as.double(Re(Ws)), i_Ds = as.double(Im(Ws)), PACKAGE = "BDgraph" )
	}
	else
	{
		result = .C( "bdmcmc_for_multi_dim", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ls), r_K = as.double(r_K), 
					i_K = as.double(i_K), as.integer(p), as.integer(Nlength), r_sigma = as.double(r_sigma), i_sigma = as.double(i_sigma), 
					r_K_hat = as.double(r_K_hat), i_K_hat = as.double(i_K_hat), p_links = as.double(p_links), as.integer(b), as.integer(b_star), 
					r_Ds = as.double(Re(Ws)), i_Ds = as.double(Im(Ws)), PACKAGE = "BDgraph" )
	}
	
	r_K_hat      = result $ r_K_hat
	i_K_hat      = result $ i_K_hat
	K_hat = matrix( complex( 1, r_K_hat, i_K_hat ), p, p * Nlength )

	last_graph = matrix( result $ G, p, p )

	last_r_K      = result $ r_K
	last_i_K      = result $ i_K
	last_K = matrix( complex( 1, last_r_K, last_i_K ), p, p * Nlength )

	if ( save.all == TRUE )
	{
		status = result $ exit
		if ( status )
		{
			mes <- paste( c( "Exit at iteration ", status ), collapse = "" )
			cat( mes, "\r" )
		}
		if ( status == 0 | status >= burnin - 1 )
		{
			size_sample_g = result $ size_sample_g
			sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
			graph_weights = result $ graph_weights[ 1 : size_sample_g ]
			all_graphs    = result $ all_graphs + 1
			all_weights   = result $ all_weights
			
			output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = K_hat, 
					all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph, last_K = last_K, status = status)
		}
		else
		{
			p_links   = matrix( result $ G, p, p ) 
			p_links[ lower.tri( p_links ) ] = 0
			output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K )
		}
	}
	else
	{
		p_links   = matrix( result $ p_links, p, p ) 
		p_links[ lower.tri( p_links ) ] = 0
		output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K )
	}

	class( output ) = "bdgraph"
	return( output )   
}
