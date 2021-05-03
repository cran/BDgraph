## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2021  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Data generator according to the graph structure                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bdgraph.sim = function( p = 10, graph = "random", n = 0, type = "Gaussian", 
						prob = 0.2, size = NULL, mean = 0, class = NULL, 
						cut = 4, b = 3, D = diag(p), K = NULL, sigma = NULL, 
						q = exp(-1), beta = 1, vis = FALSE )
{
    if( p < 2 )                       stop( " 'p' must be more than 1." )
    if( ( prob < 0 ) | ( prob > 1 ) ) stop( " Value of 'prob' must be between ( 0, 1 )." )
    if( cut < 2 )                     stop( " Value of 'cut' must be more than 1." )
    if( b <= 2 )                      stop( " Value of 'b' must be more than 2." )
    
    if( inherits( graph, "graph" ) ) graph = unclass( graph )
    
    if( is.matrix( graph ) & is.matrix( K ) ) if( nrow( graph ) != nrow( K ) ) stop( " Matrices 'graph' and 'K' have non-conforming size." )
        
    if( !is.null( size ) ) if( ( size < 0 ) | ( size > ( p * ( p - 1 ) / 2 ) ) ) stop( " Value of 'size' must be between ( 0, p*(p-1)/2 )." )
    
    if( is.matrix( K ) )
    {
        if( !isSymmetric( K ) ) stop( " Matrix 'K' must be positive definite." )
        
        graph <- "fixed"
        p     <- nrow( K )
    }
    
    if( type == "normal"     ) type = "Gaussian"
    if( type == "non-normal" ) type = "non-Gaussian"
    if( type == "discrete"   ) type = "count"
    
    if( ( type == "count" ) & ( cut == 2 ) ) type = "binary"
    
    if( is.matrix( graph ) )
	{
        if( !isSymmetric( graph ) ) stop( "Matrix 'graph' must be symmetric" )

        p = nrow( graph )
        if( ( sum( graph == 1 ) + sum( graph == 0 ) ) != ( p * p ) ) stop( "Elements of 'graph', as a matrix, must be 0 or 1." )
        
        G     <- graph
	    graph <- "fixed"
    } 
	
    # - - build the graph structure - - - - - - - - - - - - - - - - - - - - - -|
    if( sum( graph == c( "fixed", "AR1", "AR2", "circle" ) ) == 0 )
		G <- BDgraph::graph.sim( p = p, graph = graph, prob = prob, size = size, class = class )

    if( graph == "AR1" )
    {
        sigma = matrix( 0, p, p )
        
        for( i in 1 : ( p - 1 ) )
            for( j in ( i + 1 ) : p )
                sigma[ i, j ] = ( 0.7 ) ^ abs( i - j )
            
            sigma = sigma + t( sigma ) + diag( p )
            K     = solve( sigma )
            G     = 1 * ( abs( K ) > 0.02 ) 
    }
    
    if( graph == "AR2" )
    {
        K = stats::toeplitz( c( 1, 0.5, 0.25, rep( 0, p - 3 ) ) )
        G = 1 * ( abs( K ) > 0.02 ) 
    }
    
    if( graph == "circle" )
    {
        K         <- stats::toeplitz( c( 1, 0.5, rep( 0, p - 2 ) ) )
        K[ 1, p ] <- 0.4
        K[ p, 1 ] <- 0.4
        
        G = 1 * ( abs( K ) > 0.02 ) 
    }
    
    # - - generate multivariate data according to the graph structure - - - - -|
	if( n != 0 )
	{
		if( !is.null( sigma ) ) K <- solve( sigma )   

		if( is.matrix( K ) )
		{ 
			G         = 1 * ( abs( K ) > 0.02 )
			diag( G ) = 0
			# if( is.null( sigma ) ) sigma <- solve( K )	
			if( is.null( sigma ) ) sigma = stats::cov2cor( solve( K ) )
		}else{ 
			# - - Generate precision matrix according to the graph structure - |
		    if( !isSymmetric( D ) ) stop( " Matrix 'D' must be positive definite." )
		   
		    Ti        = chol( solve( D ) )
			diag( G ) = 0
			K         = matrix( 0, p, p )
			
			threshold = 1e-8  # for "rgwish_c" function in C++
	
			result = .C( "rgwish_c", as.integer(G), as.double(Ti), K = as.double(K), 
						 as.integer(b), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
			
			K     = matrix( result $ K, p, p ) 		
			# sigma = solve( K )
			sigma = stats::cov2cor( solve( K ) )
		}
		
		# - - generate multivariate normal data - - - - - - - - - - - - - - - -|
		d <- BDgraph::rmvnorm( n = n, mean = mean, sigma = sigma )
		
		# - - generate multivariate mixed data - - - - - - - - - - - - - - - - |
		not.cont = numeric( p )

		if( type == "mixed" )
		{
			# generating mixed data which are 'count', 'ordinal', 'non-Gaussian', 
			# 'binary', and 'Gaussian', respectively.
			ps = floor( p / 5 )
			
			# generating count data
			col_number        <- c( 1:ps )
			prob              <- stats::pnorm( d[, col_number] )
			d[ , col_number ] <- stats::qpois( p = prob, lambda = 10 )
			
			not.cont[ 1:ps ] = 1

			# generating ordinal data
			col_number        <- c( ( ps + 1 ):( 2 * ps ) )
			prob              <- stats::pnorm( d[ , col_number ] )
			d[ , col_number ] <- stats::qpois( p = prob, lambda = 2 )
			
			not.cont[ c( ( ps + 1 ):( 2 * ps ) ) ] = 1
			
			# generating non-Guassian data
			col_number     <- c( ( 2 * ps + 1 ):( 3 * ps ) )
			prob           <- stats::pnorm( d[ , col_number ] )
			d[,col_number] <- stats::qexp( p = prob, rate = 10 )

			# for binary data
			col_number        <- c( ( 3 * ps + 1 ):( 4 * ps ) )
			prob              <- stats::pnorm( d[ , col_number ] )
			d[ , col_number ] <- stats::qbinom( p = prob, size = 1, prob = 0.5 )
			
			not.cont[ c( ( 3 * ps + 1 ):( 4 * ps ) ) ] = 1
		}

		# - - generate multivariate continuous non-Gaussian data - - - - - - - |
		if( type == "non-Gaussian" )
		{
			# generating multivariate continuous non-Gaussian data  
			prob <- stats::pnorm( d )
			d    <- stats::qexp( p = prob, rate = 10 )
		}

		# - - generate multivariate count data - - - - - - - - - - - - - - - - |
		if( type == "count" )
		{
		    not.cont[ 1:p ] = 1
		    
			runif_m   <- matrix( stats::runif( cut * p ), nrow = p, ncol = cut )   
			marginals <- apply( runif_m, 1, function( x ) { stats::qnorm( cumsum( x / sum( x ) )[ -length( x ) ] ) } )
			if( cut == 2 ) marginals = matrix( marginals, nrow = 1, ncol = p )
				 
			for( j in 1:p )
			{
				breaks   <- c( min( d[ , j ] ) - 1, marginals[ , j ], max( d[ , j ] ) + 1 )  
				d[ , j ] <- as.integer( cut( d[ , j ], breaks = breaks, right = FALSE ) )
			}	
			
			d = d - 1
		}

		if( type == "binary" )
		{
		    not.cont[ 1:p ] = 1
		    
		    if( p > 16 ) stop( "For type 'binary', number of nodes (p) must be less than 16." )
			
			## Generate clique factors
			clique_factors = generate_clique_factors( ug = G )

			d = sample_ug( n = n, ug = G, clique_factors = clique_factors )
			
			d = d - 1
		}
		
		if( type == "dw" )
		{
            if( ( is.matrix( q ) == FALSE )  &  ( is.matrix( beta ) == FALSE) ) #q and beta can be a vector or a number
            {  
                if( length( q    ) == 1 ) q    = rep( q, time = p ) 
                if( length( beta ) == 1 ) beta = rep( beta, time = p )
                
                if( length( q    ) != p ) stop( " Length of vector 'q' has non-conforming size with 'p'. "    )
                if( length( beta ) != p ) stop( " Length of vector 'beta' has non-conforming size with 'p'. " )
                
                not.cont[ 1:p ] = 1
                Y_data <- matrix( c( 0, 1 ), nrow = n, ncol = p )
                
                while( any( apply( Y_data, 2, function( x ) { all( x %in% 0:1 ) } ) ) == TRUE ) ##detect binary variables 
                {  
                    #d <- BDgraph::rmvnorm( n = n, mean = mean, sigma = sigma )
                    d <- tmvtnorm::rtmvnorm( n = n, mean = rep( mean, p ), 
                                             sigma = sigma, lower = rep( -5, length = p ), 
                                             upper = rep( 5, length = p ) )
                    
                    
                    # STEP 2: Generate Continuous Weibull with (q, beta) - - - - - - - -|
                    #shape.c <- beta
                    #scale.c <- exp( -log( -log( q ) ) / beta )
                    
                    #for( j in 1 : p ) Y_data[ , j ] <- stats::qweibull( stats::pnorm( d[ , j ] ), shape = shape.c[ j ] , scale = scale.c[ j ] )
                    for( j in 1 : p ) 
                        Y_data[ , j ] = BDgraph::qdweibull( stats::pnorm( d[ , j ] ),  q = q[ j ], beta = beta[ j ] , zero = TRUE )		    
                    #d = Y
                    
                    # STEP 3: Transform from Continuous Weibull to Discrete Weibull - -|
                    #Y_data = floor( Y_data )
                    cat( "any binary variable", any( apply( Y_data, 2, function( x ) { all( x %in% 0 : 1 ) } ) ), "\n" )
                    #if(any( apply(Y_data, 2, function(x) { all(x %in% 0:1) })) == TRUE) tmvtnorm::rtmvnorm( n = n, mean = rep(mean, p), sigma = sigma, lower = rep(-5, length = p), upper = rep(5, length = p))
                }
            
                d = Y_data
            }
        
            if( ( is.matrix( q ) == TRUE )  &  ( is.matrix( beta ) == TRUE ) ) # q and beta come from regression 
            {
                d <- matrix( 0, nrow = n, ncol = p )
                Z <- tmvtnorm::rtmvnorm( n = n, mean = rep( mean, p ), 
                                         sigma = sigma, lower = rep( -3, length = p ), 
                                         upper = rep( 3, length = p ) )
                pnorm_Z <- stats::pnorm( Z )
                for( j in 1:p ) 
                    d[  ,j ] <- BDgraph::qdweibull( pnorm_Z[ ,j ], q = q[ , j ], beta = beta[ , j ] , zero = TRUE )
            }  
        }
		
		if( type == "NB" )
		{
            not.cont[ 1:p ] = 1
            Y.star <- matrix( c( 0, 1 ), nrow = n, ncol = p )
            while ( any( apply( Y.star, 2, function( x ) { all( x %in% 0:1 ) } ) ) == TRUE ) ##detect binary variables 
            {
                d <- tmvtnorm::rtmvnorm( n = n, mean = rep( mean, p ), sigma = sigma, 
                                         lower = rep( -5, length = p ), upper = rep( 5, length = p ) )
                
                mu   <- stats::runif( n = p, min = 0.5, max = 5 )
                size <- stats::runif( n = p, min = 0.1, max = 3 )
                
                for(j in 1 : p ) 
                    Y.star[ , j ] <- stats::qnbinom( stats::pnorm( d[ , j ] ), size = size[ j ], 
                                                     mu = mu[ j ], lower.tail = TRUE, log.p = FALSE )
                # cat( "any binary variables ", any( apply( Y , 2, function( x ) { all( x %in% 0:1 ) }) ) , "\n" )
            }
            
            d = Y.star
		}
	}
	
	# - - Saving the result - - - - - - - - - - - - - - - - - - - - - - - - - -|
	if( n != 0 )
	{
	    if( type != "dw" ){
		    simulation <- list( G = G, graph = graph, data = d, sigma = sigma, K = K, type = type, not.cont = not.cont )
	    }else{
	        simulation <- list( G = G, graph = graph, data = d, sigma = sigma, K = K, type = type, not.cont = not.cont,
	                            shape.d = beta, scale.d = q )
	    }
	}else{
		simulation <- list( G = G, graph = graph )		
	}

    # - - graph visualization - - - - - - - - - - - - - - - - - - - - - - - - -|
    if( vis == TRUE )
        BDgraph::plot.graph( G, main = "Graph structure" )
    
	class( simulation ) <- "sim"
	return( simulation )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Print function for simulation data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
print.sim = function( x, ... )
{
	p <- ncol( x $ G )
	
	if( is.null( x $ type ) )
	{
		cat( paste( "  graph generated by bdgraph.sim"                                 ), fill = TRUE )
		cat( paste( "  Graph type      =", x $ graph                                   ), fill = TRUE )
		cat( paste( "  Number of nodes =", p                                           ), fill = TRUE )
		cat( paste( "  Graph size      =", sum( x $ G ) / 2                            ), fill = TRUE )
		cat( paste( "  Sparsity        =", round( sum( x $ G ) / ( p * ( p - 1 ) ), 4) ), fill = TRUE )		
	}else{
		cat( paste( "  Data generated by bdgraph.sim"                                  ), fill = TRUE )
		cat( paste( "  Data type       =", x $ type                                    ), fill = TRUE )
		cat( paste( "  Sample size     =", nrow( x $ data )                            ), fill = TRUE )
		cat( paste( "  Graph type      =", x $ graph                                   ), fill = TRUE )
		cat( paste( "  Number of nodes =", p                                           ), fill = TRUE )
		cat( paste( "  Graph size      =", sum( x $ G ) / 2                            ), fill = TRUE )
		cat( paste( "  Sparsity        =", round( sum( x $ G ) / ( p * ( p - 1 ) ), 4) ), fill = TRUE )
	}
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# plot for class "sim" from bdgraph.sim function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
plot.sim = function( x, ... )
{
    BDgraph::plot.graph( x, ... )
}		
       
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Function for exact sampling from binary data
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
sample_ug = function( n = 1, ug = diag( 3 ), clique_factors = NULL )
{
	p = ncol( ug ) # p smaller than 17 check
	if( p > 16 ) stop( "number of nodes must be smaller than 16." )
	ug[ lower.tri( ug, diag = TRUE ) ] = 0
	if( is.null( clique_factors ) ) clique_factors = generate_clique_factors( ug )
	
	prob = calc_joint_dist( ug, clique_factors )
	noc  = length( prob )
	ind  = sample( 1:noc, n, replace = TRUE, prob = prob )
	oc   = sapply( 1:noc, function( x ){ as.integer( intToBits( x ) ) } )
	oc   = 2 - t( oc[ 1:p, ] )
	data = oc[ ind, ]
	
	return( data )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
generate_clique_factors = function( ug )
{
	ug[ lower.tri( ug, diag = TRUE ) ] = 0   
	p              = ncol( ug )
	edges          = which( ug == 1, arr.ind = T )
	a              = igraph::make_undirected_graph( c( t( edges ) ), p )
	cliques        = igraph::max_cliques( a )
	clique_factors = vector( 'list', length( cliques ) )
	
	for ( i in 1:length( cliques ) )
	{
		clq                 = cliques[[i]]
		clique_factors[[i]] = stats::runif( 2 ^ length( clq ) )
	}
	
	return( clique_factors )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
calc_joint_dist = function( ug, clique_factors )
{
	p          = ncol( ug )
	oc         = sapply( 1:( 2 ^ p ), function( x ){ as.integer( intToBits( x ) ) } )
	oc         = 2 - t( oc[ 1:p, ] )

	edges      = which( ug == 1, arr.ind = T )
	a          = igraph::make_undirected_graph( c( t( edges ) ), p )

	joint_dist = rep( 1, 2 ^ p )
	cliques    = igraph::max_cliques( a )
	
	for ( i in 1:length( cliques ) )
	{
		clq        = cliques[[i]]
		k          = length( clq )
		temp       = sapply( 1:( 2 ^ k ), function( x ){ as.integer( intToBits( x ) ) } )
		clq_oc     = 2 - t( temp[1:k, ] )
		clq_factor = clique_factors[[i]]
		
		for ( j in 1:nrow( clq_oc ) )
		{
			oc_col_clq = oc[ , clq ]
			if( is.null( dim( oc_col_clq ) ) ) 
				oc_col_clq = matrix( oc_col_clq, nrow = length( oc_col_clq ), ncol = 1 )
		
			ind             = apply( oc_col_clq, 1, identical, clq_oc[ j, ] )
			joint_dist[ind] = joint_dist[ ind ] * clq_factor[ j ]
		}				
	}	
	
	joint_dist = joint_dist / sum( joint_dist )
	
	return( joint_dist )	
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |


