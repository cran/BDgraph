## ------------------------------------------------------------------------------------------------|
#     Copyright (C) 2012 - 2018  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of BDgraph package.                                                        |
#                                                                                                  |
#     BDgraph is free software: you can redistribute it and/or modify it under                     |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## ------------------------------------------------------------------------------------------------|
#     Data generator according to the graph structure                                              |
## ------------------------------------------------------------------------------------------------|

bdgraph.sim = function( p = 10, graph = "random", n = 0, type = "Gaussian", 
						prob = 0.2, size = NULL, mean = 0, class = NULL, 
						cut = 4, b = 3, D = diag(p), K = NULL, sigma = NULL, 
						vis = FALSE )
{
    if( p < 2 )                       stop( " 'p' must be more than 1" )
    if( ( prob < 0 ) | ( prob > 1 ) ) stop( " Value of 'prob' must be between ( 0, 1 )" )
    if( cut < 2 )                     stop( " Value of 'cut' must be more than 1" )
    if( b <= 2 )                      stop( " Value of 'b' must be more than 2" )
    
    if( is.matrix( graph ) & is.matrix( K ) ) if( nrow( graph ) != nrow( K ) ) stop( " matrices 'graph' and 'K' have non-conforming size" )
        
    if( !is.null( size ) ) if( ( size < 0 ) | ( size > ( p * ( p - 1 ) / 2 ) ) ) stop( " Value of 'size' must be between ( 0, p*(p-1)/2 )" )
    
    if( is.matrix( K ) )
    {
        if( !isSymmetric( K ) ) stop( " Matrix 'K' must be positive definite" )
        
        graph <- "fixed"
        p     <- nrow( K )
    }
    
    if( ( type == "discrete" ) & ( cut == 2 ) ) type = "binary"
    
    if( type == "normal"     ) type = "Gaussian"
    if( type == "non-normal" ) type = "non-Gaussian"
    
    if( is.matrix( graph ) )
	{
        if( !isSymmetric( graph ) ) stop( "Matrix 'graph' must be symmetric" )

        p = nrow( graph )
        if( ( sum( graph == 1 ) + sum( graph == 0 ) ) != ( p * p ) ) stop( "Element of 'graph', as a matrix, must be 0 or 1." )
        
        G     <- graph
	    graph <- "fixed"
    } 
	
    #--- build the graph structure ----------------------------------------------------------------|
    if( sum( graph == c( "fixed", "AR1", "AR2", "star" ) ) == 0 )
		G <- BDgraph::graph.sim( p = p, graph = graph, prob = prob, size = size, class = class )

    if( graph == "AR1" )
    {
        sigma = matrix( 0, p, p )
        
        for( i in 1 : ( p - 1 ) )
            for( j in ( i + 1 ) : p )
                sigma[i, j] = ( 0.7 ) ^ abs( i - j )
            
            sigma = sigma + t( sigma ) + diag( p )
            K     = solve( sigma )
            G     = 1 * ( abs( K ) > 0.02 ) 
    }
    
    if( graph == "AR2" )
    {
        K = stats::toeplitz( c( 1, 0.5, 0.25, rep( 0, p - 3 ) ) )
        G = 1 * ( abs( K ) > 0.02 ) 
    }
    
    if( graph == "star" )
    {
        K                 <- diag( p )
        K[ 1, ( 2 : p ) ] <- 0.1
        K[ ( 2 : p ), 1 ] <- 0.1
        
        G = 1 * ( abs( K ) > 0.02 ) 
    }

	#--- generate multivariate data according to the graph structure ------------------------------|
	if( n != 0 )
	{
		if( !is.null( sigma ) ) K <- solve( sigma )   

		if( is.matrix( K ) )
		{ 
			G         <- 1 * ( abs( K ) > 0.02 )
			diag( G ) <- 0
			if( is.null( sigma ) ) sigma <- solve( K )	
		}else{ 
			#--- Generate precision matrix according to the graph structure -----------------------|
		    if( !isSymmetric( D ) ) stop( " Matrix 'D' must be positive definite" )
		   
		    Ti        = chol( solve( D ) )
			diag( G ) = 0
			K         = matrix( 0, p, p )
			
			threshold = 1e-8  # for "rgwish_c" function in C++
	
			result = .C( "rgwish_c", as.integer(G), as.double(Ti), K = as.double(K), 
						 as.integer(b), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
			
			K     = matrix( result $ K, p, p ) 		
			sigma = solve( K )
		}
		
		#--- generate multivariate normal data ----------------------------------------------------|
		d <- BDgraph::rmvnorm( n = n, mean = mean, sigma = sigma )
		
		#--- generate multivariate mixed data -----------------------------------------------------|
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

		#--- generate multivariate continuous non-Gaussian data -----------------------------------|
		if( type == "non-Gaussian" )
		{
			# generating multivariate continuous non-Gaussian data  
			prob <- stats::pnorm( d )
			d    <- stats::qexp( p = prob, rate = 10 )
		}

		#--- generate multivariate discrete data --------------------------------------------------|
		if( type == "discrete" )
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
		    
		    if( p > 16 ) stop( "For type 'binary', number of nodes (p) must be less than 16" )
			
			## Generate clique factors
			clique_factors = generate_clique_factors( ug = G )

			d = sample_ug( n = n, ug = G, clique_factors = clique_factors )
			
			d = d - 1
		}
	}
	
	#--- Saving the result ------------------------------------------------------------------------|
	if( n != 0 )
	{
		simulation <- list( G = G, graph = graph, data = d, sigma = sigma, K = K, type = type, not.cont = not.cont )
	}else{
		simulation <- list( G = G, graph = graph )		
	}

    #--- graph visualization ----------------------------------------------------------------------|
    if( vis == TRUE )
    {
        true_graph = as.matrix( G )
        graphG <- igraph::graph.adjacency( true_graph, mode = "undirected", diag = FALSE )
        
        if( p < 20 ) size = 10 else size = 2
        igraph::plot.igraph( graphG, layout = igraph::layout.circle, main = "Graph structure", 
                     vertex.color = "white", vertex.size = size, vertex.label.color = 'black' )
    }
    
	class( simulation ) <- "sim"
	return( simulation )
}
    
## ------------------------------------------------------------------------------------------------|
# Print function for simulation data
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
    
## ------------------------------------------------------------------------------------------------|
# plot for class "sim" from bdgraph.sim function
plot.sim = function( x, main = NULL, layout = layout.circle, ... )
{
    true_graph = as.matrix( x $ G )
    if( is.null( main ) ) main = "Graph structure"
  	g_igraph <- igraph::graph.adjacency( true_graph, mode = "undirected", diag = FALSE )
	
  	igraph::plot.igraph( g_igraph, main = main, layout = layout, ... )
}		
       
## ------------------------------------------------------------------------------------------------|
# Function for exact sampling from binary data
## ------------------------------------------------------------------------------------------------|
sample_ug = function( n = 1, ug = diag( 3 ), clique_factors = NULL )
{
	p = ncol( ug ) # p smaller than 17 check
	if( p > 16 ) stop( "number of nodes must be smaller than 16" )
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
    
## ------------------------------------------------------------------------------------------------|
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
    
## ------------------------------------------------------------------------------------------------|
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
    


