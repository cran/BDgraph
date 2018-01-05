## ----------------------------------------------------------------------------|
# Main function of BDgraph package: BDMCMC algorithm for graphical models 
## ----------------------------------------------------------------------------|
bdgraph = function( data, n = NULL, method = "ggm", algorithm = "bdmcmc", 
					iter = 5000, burnin = iter / 2, g.start = "empty", g.space = NULL, g.prior = 0.5, 
					prior.df = 3, multi.update = NULL, save.all = FALSE, print = 1000, cores = "all" )
{
	check.os( os = 2 )	
	
	if( cores == "all" ) cores = detect_cores()
	
	tmp   <- .C( "check_nthread", cores = as.integer(cores), PACKAGE = "BDgraph" )
	cores <- tmp $ cores
	
	.C( "omp_set_num_cores", as.integer( cores ), PACKAGE = "BDgraph" )
	
	burnin <- floor( burnin )
	
	if( class( data ) == "sim" ) data <- data $ data

	if( !is.matrix( data ) & !is.data.frame( data ) ) stop( "Data should be a matrix or dataframe" )
	if( is.data.frame( data ) ) data <- data.matrix( data )
	if( iter < burnin )   stop( "Number of iteration must be more than number of burn-in" )

	if( any( is.na( data ) ) ) 
	{
		if( method == "ggm" ) { stop( "ggm method does not deal with missing values. You could choose option method = gcgm" ) }	
		gcgm_NA = 1
	}else{
		gcgm_NA = 0
	}
		
	dimd <- dim( data )
	p    <- dimd[2]
	if( is.null( n ) ) n <- dimd[1]

	if( !is.matrix( g.prior ) )
	{
		if( ( g.prior <= 0 ) | ( g.prior >= 1 ) ){
			stop( "g.prior must be between 0 and 1" )
		}else{
			g.prior = matrix( g.prior, p, p )
		}
	}
	g_prior = g.prior

	if( method == "gcgm" )
	{
		if( isSymmetric( data ) ) stop( "method='gcgm' requires all data" )
		
 		R <- 0 * data
		for( j in 1:p ) R[,j] = match( data[ , j], sort( unique( data[ , j] ) ) ) 
		R[ is.na(R) ] = 0     # dealing with missing values	

		# copula for continuous non-Gaussian data
		if( gcgm_NA == 0 && min( apply( R, 2, max ) ) > ( n - 5 * n / 100 ) )
		{
			# copula transfer 
			data = qnorm( apply( data, 2, rank ) / ( n + 1 ) )
#~ 			data = data / sd( data[ , 1] )
		
			method = "ggm"
		}else{	# for non-Gaussian data
			Z              <- qnorm( apply( data, 2, rank, ties.method = "random" ) / ( n + 1 ) )
			Zfill          <- matrix( rnorm( n * p ), n, p )   # for missing values
			Z[is.na(data)] <- Zfill[is.na(data)]               # for missing values
			Z              <- t( ( t(Z) - apply( Z, 2, mean ) ) / apply( Z, 2, sd ) )
			S              <- t(Z) %*% Z
		}
	} 
	
	if( method == "ggm" ) 
	{
		if( isSymmetric( data ) )
		{
			if ( is.null(n) ) stop( "Please specify the number of observations 'n'" )
			cat( "Input is identified as the covriance matrix. \n" )
			S <- data
		}else{
 			S <- t(data) %*% data
		}
	}
   
	D      = diag( p )
	b      = prior.df
	b_star = b + n
	Ds     = D + S
	Ts     = chol( solve( Ds ) )
	Ti     = chol( solve( D ) )   # only for double Metropolic-Hastings algorithms 

	if( class( g.start ) == "bdgraph" ) 
	{
		G <- g.start $ last_graph
		K <- g.start $ last_K
	} 

	if( class( g.start ) == "sim" ) 
	{
		G <- as.matrix( g.start $ G )
		K <- as.matrix( g.start $ K )
	} 
	
	if( class( g.start ) == "character" && g.start == "empty"  )
	{
		G = matrix( 0, p, p )
		K = G
		
		result = .C( "rgwish_c", as.integer(G), as.double(Ts), K = as.double(K), as.integer(b_star), as.integer(p), PACKAGE = "BDgraph" )
		K      = matrix ( result $ K, p, p ) 
	}
	
	if( class( g.start ) == "character" && g.start == "full" )
	{
		G         = matrix( 1, p, p )
		diag( G ) = 0
		K         = matrix( 0, p, p )

		result = .C( "rwish_c", as.double(Ts), K = as.double(K), as.integer(b_star), as.integer(p), PACKAGE = "BDgraph" )
		K      = matrix ( result $ K, p, p ) 
	}	

	if( is.matrix( g.start ) )
	{
		G       = g.start
		diag(G) = 0
		
		K      = matrix( 0, p, p )
		result = .C( "rgwish_c", as.integer(G), as.double(Ts), K = as.double(K), as.integer(b_star), as.integer(p), PACKAGE = "BDgraph" )
		K      = matrix ( result $ K, p, p ) 	
	}
			
	if( save.all == TRUE )
	{
		qp1           = ( p * ( p - 1 ) / 2 ) + 1
		string_g      = paste( c( rep( 0, qp1 ) ), collapse = '' )
		sample_graphs = c( rep ( string_g, iter - burnin ) )  # vector of numbers like "10100" 
		graph_weights = c( rep ( 0, iter - burnin ) )         # waiting time for every state
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
	
	K_hat      = matrix( 0, p, p )
	last_graph = K_hat
	last_K     = K_hat

	if( ( is.null( multi.update ) ) && ( p > 10 & iter > ( 5000 / p ) ) )
		multi.update = floor( p / 10 )
	
	if( is.null( multi.update ) ) multi.update = 1
	multi_update = multi.update

	if( ( p < 10 ) && ( multi_update > 1 ) )      cat( " WARNING: for the cases p < 10, multi.update must be 1 " )
	if( multi_update > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: multi.update must be smaller " )
	
	mes <- paste( c( iter, " iteration is started.                    " ), collapse = "" )
	cat( mes, "\r" )

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
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "ggm_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "ggm_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(b), as.integer(b_star), as.double(Ds), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "ggm_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 ) )
		{
			result = .C( "gcgm_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "gcgm_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}
   
		if( ( method == "gcgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "gcgm_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}	

		# for Double Metropolis-Hasting 
		if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( multi_update == 1 ) )
		{
			result = .C( "ggm_DMH_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "ggm_DMH_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "rj-dmh" ) )
		{
			result = .C( "ggm_DMH_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), 
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( multi_update == 1 ) )
		{
			result = .C( "gcgm_DMH_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}
 
 		if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( multi_update != 1 ) )
		{
			counter_all_g = 0
			result = .C( "gcgm_DMH_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "rj-dmh" ) )
		{
			result = .C( "gcgm_DMH_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
						sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}	
      
	}else{
		
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 )  )
		{
			result = .C( "ggm_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), 
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			result = .C( "ggm_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), 
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(Ds), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}		
    
		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "ggm_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), 
						K_hat = as.double(K_hat), p_links = as.integer(p_links),
						as.integer(b), as.integer(b_star), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update == 1 )  )
		{
			result = .C( "gcgm_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}
    
		if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( multi_update != 1 ) )
		{
			result = .C( "gcgm_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "gcgm_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						K_hat = as.double(K_hat), p_links = as.integer(p_links),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}	

		# for Double Metropolis-Hasting 
		if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( multi_update == 1 )  )
		{
			result = .C( "ggm_DMH_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), 
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( multi_update != 1 ) )
		{
			result = .C( "ggm_DMH_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), 
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}		
    
		if( ( method == "ggm" ) && ( algorithm == "rj-dmh" ) )
		{
			result = .C( "ggm_DMH_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), 
						K_hat = as.double(K_hat), p_links = as.integer(p_links),
						as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(print), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( multi_update == 1 )  )
		{
			result = .C( "gcgm_DMH_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( multi_update != 1 ) )
		{
			result = .C( "gcgm_DMH_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						K_hat = as.double(K_hat), p_links = as.double(p_links),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(multi_update), as.integer(print), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "rj-dmh" ) )
		{
			result = .C( "gcgm_DMH_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.integer(g_space), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						K_hat = as.double(K_hat), p_links = as.integer(p_links),
						as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(print), PACKAGE = "BDgraph" )
		}	
	}
## ----------------------------------------------------------------------------|

	label      = colnames( data )

	K_hat      = matrix( result $ K_hat, p, p, dimnames = list( label, label ) ) 
	last_graph = matrix( result $ G    , p, p, dimnames = list( label, label ) )
	last_K     = matrix( result $ K    , p, p )

#	colnames( last_graph ) = colnames( data )

	if( save.all == TRUE )
	{
		if( algorithm == "rjmcmc" ) K_hat = K_hat / ( iter - burnin )		
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

		output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = K_hat, 
					all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph, last_K = last_K )
	}else{
#		label  = colnames( last_graph )
		p_links = matrix( result $ p_links, p, p, dimnames = list( label, label ) ) 

		if( ( algorithm == "rjmcmc" ) | ( algorithm == "rj-dmh" ) )
		{
			p_links = p_links / ( iter - burnin )
			K_hat   = K_hat / ( iter - burnin )
		}
		p_links[ lower.tri( p_links ) ] = 0
#		colnames( p_links ) = colnames( data )
		output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K )
	}
	
	class( output ) = "bdgraph"
	return( output )   
}
      
## ----------------------------------------------------------------------------|
# summary of bdgraph output
summary.bdgraph = function( object, round = 2, vis = TRUE, ... )
{
	p_links    = object $ p_links
	p          = nrow( object $ last_graph )
	label      = colnames( object $ last_graph )
	selected_g = matrix( 0, p, p, dimnames = list( label, label ) )	

	if( !is.null( object $ graph_weights ) )
	{
		sample_graphs = object $ sample_graphs
		graph_weights = object $ graph_weights
		max_gWeights  = max( graph_weights )
		sum_gWeights  = sum( graph_weights )
		max_prob_G    = max_gWeights / sum_gWeights

		if ( is.null( label ) ) label <- as.character( 1 : p )
		vec_G    <- c( rep( 0, p * ( p - 1 ) / 2 ) )		
		indG_max <- sample_graphs[ which( graph_weights == max_gWeights ) ]
		vec_G[ which( unlist( strsplit( as.character( indG_max ), "" ) ) == 1 ) ] = 1
		selected_g[ upper.tri( selected_g ) ] <- vec_G 
	}else{
		selected_g[ p_links > 0.5 ]  = 1
		selected_g[ p_links <= 0.5 ] = 0
	}

	if( vis )
	{
		# plot selected graph (graph with the highest posterior probability)
		G  <- graph.adjacency( selected_g, mode = "undirected", diag = FALSE )
		 
		if( !is.null( object $ graph_weights ) ) 
		{
			op       = par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 
			subGraph = paste( c( "Posterior probability = ", max_prob_G ), collapse = "" )
		}else{
			subGraph = "Selected graph with edge posterior probability = 0.5"
		}
			
		if( p < 20 ) size = 15 else size = 2
		plot.igraph( G, layout = layout.circle, main = "Selected graph", sub = subGraph, vertex.color = "white", vertex.size = size, vertex.label.color = 'black' )
		
		if( !is.null( object $ graph_weights ) )
		{
			# plot posterior distribution of graph
			plot( x = 1 : length( graph_weights ), y = graph_weights / sum_gWeights, type = "h", main = "Posterior probability of graphs",
				 ylab = "Pr(graph|data)", xlab = "graph" )
			
			abline( h = max_prob_G, col = "red" )
			text( which( max_gWeights == graph_weights )[1], max_prob_G, "Pr(selected graph|data)", col = "gray60", adj = c( 0, +1 ) )
			
			# plot posterior distribution of graph size
			sizesample_graphs = sapply( sample_graphs, function( x ) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )
			xx       <- unique( sizesample_graphs )
			weightsg <- vector()

			for( i in 1 : length(xx) ) weightsg[i] <- sum( graph_weights[ which( sizesample_graphs == xx[i] ) ] )

			plot( x = xx, y = weightsg / sum_gWeights, type = "h", main = "Posterior probability of graphs size", ylab = "Pr(graph size|data)", xlab = "Graph size" )

			# plot trace of graph size
			all_graphs     = object $ all_graphs
			sizeall_graphs = sizesample_graphs[ all_graphs ]
			  
			plot( x = 1 : length( all_graphs ), sizeall_graphs, type = "l", main = "Trace of graph size", ylab = "Graph size", xlab = "Iteration" )
			
			abline( h = sum( selected_g ), col = "red" )	  
			
			par( op )
		}
	}
	
	# p_links
	if( !is.null( object $ graph_weights ) )
	{
		pvec <- 0 * vec_G
		for( i in 1 : length( sample_graphs ) )
		{
			which_edge       <- which( unlist( strsplit( as.character( sample_graphs[i] ), "" ) ) == 1 )
			pvec[which_edge] <- pvec[which_edge] + graph_weights[i]
		}
		p_links                     <- 0 * selected_g
		p_links[upper.tri(p_links)] <- pvec / sum_gWeights
	}
	
	K_hat = object $ K_hat
	
	if( is.null( K_hat ) )			  
		return( list( selected_g = Matrix( selected_g, sparse = TRUE ), p_links = Matrix( round( p_links, round ), sparse = TRUE ) ) )
	else
		return( list( selected_g = Matrix( selected_g, sparse = TRUE ), p_links = Matrix( round( p_links, round ), sparse = TRUE ), K_hat = round( K_hat, round ) ) )
}  
   
## ----------------------------------------------------------------------------|
# plot for class bdgraph
plot.bdgraph = function( x, cut = NULL, number.g = 1, layout = layout.circle, ... )
{
	if( !is.matrix( x ) )
	{
		if( is.null( cut ) ) cut = 0.5
		
		if( is.null( cut ) )
		{
			sample_graphs = x $ sample_graphs
			graph_weights = x $ graph_weights
			prob_G        = graph_weights / sum( graph_weights )
			sort_prob_G   = sort( prob_G, decreasing = TRUE )

			p             = nrow( x $ last_graph )
			label         = colnames( x $ last_graph )
			list_G        = replicate( number.g, matrix( 0, p, p, dimnames = list( label, label ) ), simplify = FALSE )
			vec_G         = c( rep( 0, p * ( p - 1 ) / 2 ) )

			if( number.g == 2 ) op <- par( mfrow = c( 1, 2 ), pty = "s" )
			if( number.g > 2 & number.g < 7 )  op <- par( mfrow = c( 2, number.g %% 2 + trunc( number.g / 2 ) ), pty = "s" )

			for( i in 1 : number.g )
			{
				if( number.g > 6 ) dev.new()  
				indG_i <- sample_graphs[ which( prob_G == sort_prob_G[i] )[1] ]
				vec_G  <- 0 * vec_G
				vec_G[ which( unlist( strsplit( as.character(indG_i), "" ) ) == 1 ) ] <- 1
				list_G[[i]][ upper.tri( list_G[[i]] ) ] <- vec_G
				G    <- graph.adjacency( list_G[[i]], mode = "undirected", diag = FALSE )
				main <- ifelse( i == 1, "Graph with highest probability", paste( c( i, "th graph" ), collapse = "" ) )
				plot.igraph( G, layout = layout, main = main, sub = paste( c( "Posterior probability = ", 
							round( sort_prob_G[i], 6 ) ), collapse = "" ), ... )	   
			}
			
			if( number.g > 1 & number.g < 7 ) par( op )
		}else{
			if( ( cut < 0 ) || ( cut > 1 ) ) stop( " Value of 'cut' must be between 0 and 1. " )

			p_links = x $ p_links
			if( is.null( p_links ) ) p_links = plinks( x )
			selected_g                   = 0 * p_links
			selected_g[ p_links > cut ]  = 1
			selected_g[ p_links <= cut ] = 0		

			G = graph.adjacency( selected_g, mode = "undirected", diag = FALSE )
			plot.igraph( G, layout = layout, main = "Selected graph", sub = paste0( "Edge posterior probability = ", cut ), ... )	   		
		}
	}else{
			G = graph.adjacency( x, mode = "undirected", diag = FALSE )
			plot.igraph( G, layout = layout, main = "Graph with highest posterior probability", ... )	   		
	}
}
     
## ----------------------------------------------------------------------------|
# print of the bdgraph output
print.bdgraph = function( x, round = 2, ... )
{
	if( !is.matrix( x ) )
	{
		p_links = x $ p_links
		
		if( !is.null( x $ graph_weights ) )
		{
			p             = nrow( x $ last_graph )
			sample_graphs  = x $ sample_graphs
			graph_weights  = x $ graph_weights
			# selected graph
			max_gWeights  = max( graph_weights )
			sum_gWeights  = sum( graph_weights )
			vec_G         = c( rep( 0, p * ( p - 1 ) / 2 ) )
			indG_max      = sample_graphs[ which( graph_weights == max_gWeights )[1] ]
			vec_G[ which( unlist( strsplit( as.character( indG_max ), "" ) ) == 1 ) ] = 1

			label      = colnames( x $ last_graph )
			selected_g = matrix( 0, p, p, dimnames = list( label, label ) )	
			selected_g[upper.tri(selected_g)] = vec_G
		
		}else{
			selected_g                   = 0 * p_links
			selected_g[ p_links > 0.5 ]  = 1
			selected_g[ p_links <= 0.5 ] = 0	
		}
		
	}else{
		selected_g = unclass( x )
	}
	
	cat( paste( "" ), fill = TRUE )
	cat( paste( "Adjacency matrix of selected graph" ), fill = TRUE )
	cat( paste( "" ), fill = TRUE )
	
	if( !is.matrix( x ) )
	{	
		printSpMatrix( Matrix( selected_g, sparse = TRUE ), col.names = TRUE, note.dropping.colnames = FALSE )
		cat( paste( "" ), fill = TRUE )
		cat( paste( "Size of selected graph = ", sum( selected_g ) ), fill = TRUE )
	
	}else{
		print( selected_g )
		cat( paste( "" ), fill = TRUE )
		cat( paste( "Size of selected graph = ", sum( selected_g ) / 2 ), fill = TRUE )	
	}

	if( !is.matrix( x ) )
	{
		if( !is.null( x $ graph_weights ) )
			cat( paste( "Posterior probability of selected graph = ", max_gWeights / sum_gWeights ), fill = TRUE )  
		else
			cat( paste( "Edge posterior probability of selected graph = ", 0.5 ), fill = TRUE )
	}
	
	cat( paste( "" ), fill = TRUE )
} 
   




