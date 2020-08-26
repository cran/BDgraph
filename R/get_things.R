## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2020  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  Get graph from class objects "sim", "graph", "bdgraph", "ssgraph", or "select"
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_graph = function( obj_G, cut = 0.5 )
{
    if( is.matrix( obj_G ) ) 
    {
        if( nrow( obj_G ) != ncol( obj_G ) ) stop( "Adjacency matrix must be squere." )
        if( ( sum( obj_G == 0 ) + sum( obj_G == 1 ) ) != ( ncol( obj_G ) ^ 2 ) ) stop( "Elements of matrix must be 0 or 1." )
        G = unclass( obj_G )
    }else{
        if(   inherits( obj_G, "sim"     ) ) G <- unclass( obj_G $ G ) 
        if(   inherits( obj_G, "graph"   ) ) G <- unclass( obj_G ) 
        
        if( ( inherits( obj_G, "bdgraph" ) ) | ( inherits( obj_G, "ssgraph" ) ) ) G <- BDgraph::select( obj_G, cut = cut ) 
        if(   inherits( obj_G, "select"  ) ) G <- obj_G $ refit
        
        G = as.matrix( G )
    }
    
    return( G )    
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_g_prior = function( g.prior, p )
{
    if( is.data.frame( g.prior ) ) g.prior <- data.matrix( g.prior )
    if( inherits( g.prior, "dtCMatrix" ) ) g.prior = as.matrix( g.prior )
    if( ( inherits( g.prior, "bdgraph"  ) ) | ( inherits( g.prior, "ssgraph" ) ) ) g.prior <- BDgraph::plinks( g.prior )
    if( inherits( g.prior, "sim" ) ) 
    {
        K       = as.matrix( g.prior $ K )
        g.prior = abs( K / diag( K ) )
    }
    
    if( !is.matrix( g.prior ) )
    {
        if( ( g.prior <= 0 ) | ( g.prior >= 1 ) ) stop( " 'g.prior' must be between 0 and 1." )
        g.prior = matrix( g.prior, p, p )
    }else{
        if( ( nrow( g.prior ) != p ) | ( ncol( g.prior ) != p ) ) stop( " 'g.prior' and 'data' have non-conforming size." )
        if( any( g.prior < 0 ) || any( g.prior > 1 ) ) stop( " Elements of 'g.prior', as a matrix, must be between 0 and 1." )
    }
    
    g.prior[ lower.tri( g.prior, diag = TRUE ) ] <- 0
    g.prior = g.prior + t( g.prior )
    
    return( g.prior )
}
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_g_start = function( g.start, g_prior, p )
{
    if( is.matrix( g.start ) )
    {
        if( ( sum( g.start == 0 ) + sum( g.start == 1 ) ) != ( p ^ 2 ) ) stop( " Elements of 'g.start', as a matrix, must be 0 or 1." )
        G = g.start
    }
    
    if( inherits( g.start, "sim"   ) ) G <- unclass( g.start $ G )
    if( inherits( g.start, "graph" ) ) G <- unclass( g.start )
    
    if( ( inherits( g.start, "bdgraph"   ) ) |  ( inherits( g.start, "ssgraph" ) ) ) G <- g.start $ last_graph
    if( ( inherits( g.start, "character" ) ) && ( g.start == "empty" ) ) G = matrix( 0, p, p )
    if( ( inherits( g.start, "character" ) ) && ( g.start == "full"  ) ) G = matrix( 1, p, p )
    
    if( ( nrow( G ) != p ) | ( ncol( G ) != p ) ) stop( " 'g.start' and 'data' have non-conforming size." )
    
    G[ g_prior == 1 ] = 1
    G[ g_prior == 0 ] = 0
    
    G[ lower.tri( G, diag( TRUE ) ) ] <- 0
    G  = G + t( G )
    
    return( G = G )
}
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_K_start = function( G, g.start, Ts, b_star, threshold )
{
    p = ncol( G )
    
    if( ( inherits( g.start, "bdgraph" ) ) | ( inherits( g.start, "ssgraph" ) ) ) 
        K <- g.start $ last_K
    
    if( inherits( g.start, "sim" ) )  K <- g.start $ K
    
    if( ( !inherits( g.start, "bdgraph" ) ) && ( !inherits( g.start, "ssgraph" ) ) && ( !inherits( g.start, "sim" ) ) )
    {
        K = G
        
        result = .C( "rgwish_c", as.integer(G), as.double(Ts), K = as.double(K), as.integer(b_star), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
        K      = matrix ( result $ K, p, p ) 
    }
    
    return( K = K )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_S_n_p = function( data, method, n, not.cont = NULL )
{
    if( inherits( data, "sim" ) )
    {
        not.cont <- data $ not.cont  # Do not change the order of these links
        data     <- data $ data
    }
    
    if( !is.matrix( data ) & !is.data.frame( data ) ) stop( " Data must be a matrix or dataframe." )
    if( is.data.frame( data ) ) data <- data.matrix( data )
    
    if( any( is.na( data ) ) ) 
    {
        if( method == "dw"  ) stop( " 'bdgraph.dw()' does not deal with missing values." )	
        if( method == "ggm" ) stop( " 'ggm' method does not deal with missing values. You could choose option method = gcgm." )	
        gcgm_NA = 1
    }else{
        gcgm_NA = 0
    }
    
    if( isSymmetric( data ) )
    {
        if( method == "gcgm" ) stop( " method='gcgm' requires all data." )
        if( method == "dw"   ) stop( " method='dw' requires all data." )
        if( is.null( n )     ) stop( " Please specify the number of observations 'n'." )
    }
    
    p <- ncol( data )
    if( p < 3 ) stop( " Number of variables/nodes ('p') must be more than 2." )
    if( is.null( n ) ) n <- nrow( data )
    
    if( method == "ggm" ) 
    {
        if( isSymmetric( data ) )
        {
            cat( "Input is identified as the covariance matrix. \n" )
            S <- data
        }else{
            S <- t( data ) %*% data
        }
    }
    
    if( method == "gcgm" )
    {
        if( is.null( not.cont ) )
        {
            not.cont = c( rep( 1, p ) )
            for( j in 1:p )
                if( length( unique( data[ , j ] ) ) > min( n / 2 ) ) not.cont[ j ] = 0
        }else{
            if( !is.vector( not.cont )  ) stop( " 'not.cont' must be a vector with length of number of variables." )
            if( length( not.cont ) != p ) stop( " 'not.cont' must be a vector with length of number of variables." )
            if( ( sum( not.cont == 0 ) + sum( not.cont == 1 ) ) != p ) stop( " Elements of 'not.cont', as a vector, must be 0 or 1." )
        }
        
        R <- 0 * data
        for( j in 1:p )
            if( not.cont[ j ] )
                R[ , j ] = match( data[ , j ], sort( unique( data[ , j ] ) ) ) 
        R[ is.na( R ) ] = -1000     # dealing with missing values	
        
        # copula for continuous non-Gaussian data
        if( ( gcgm_NA == 0 ) && ( min( apply( R, 2, max ) ) > ( n - 5 * n / 100 ) ) )
        {
            # copula transfer 
            data = stats::qnorm( apply( data, 2, rank ) / ( n + 1 ) )
            data = t( ( t( data ) - apply( data, 2, mean ) ) / apply( data, 2, stats::sd ) )
            
            method = "ggm"
        }else{	
            # for non-Gaussian data
            Z                  <- stats::qnorm( apply( data, 2, rank, ties.method = "random" ) / ( n + 1 ) )
            Zfill              <- matrix( stats::rnorm( n * p ), n, p )   # for missing values
            Z[ is.na( data ) ] <- Zfill[ is.na( data ) ]                  # for missing values
            Z                  <- t( ( t( Z ) - apply( Z, 2, mean ) ) / apply( Z, 2, stats::sd ) )
            S                  <- t( Z ) %*% Z
        }
    } 
    
    if( method == "dw" )
    {
        if( is.null( not.cont ) )
        {
            not.cont = c( rep( 1, p ) )
        }else{
            if( !is.vector( not.cont )  ) stop( " 'not.cont' must be a vector with length of number of variables." )
            if( length( not.cont ) != p ) stop( " 'not.cont' must be a vector with length of number of variables." )
            if( ( sum( not.cont == 0 ) + sum( not.cont == 1 ) ) != p ) stop( " Elements of 'not.cont', as a vector, must be 0 or 1." )
        }
        
        # for non-Gaussian data
        Z                  <- stats::qnorm( apply( data, 2, rank, ties.method = "random" ) / ( n + 1 ) )
        Zfill              <- matrix( stats::rnorm( n * p ), n, p )   # for missing values
        Z[ is.na( data ) ] <- Zfill[ is.na( data ) ]                  # for missing values
        Z                  <- t( ( t( Z ) - apply( Z, 2, mean ) ) / apply( Z, 2, stats::sd ) )
        S                  <- t( Z ) %*% Z
    } 
    
    if( method == "ggm" ) 
        list_out = list( method = method, S = S, n = n, p = p, colnames_data = colnames( data ) )

    if( method == "gcgm" ) 
        list_out = list( method = method, S = S, n = n, p = p, colnames_data = colnames( data ), not.cont = not.cont, R = R, Z = Z, data = data, gcgm_NA = gcgm_NA )
    
    if( method == "dw" ) 
        list_out = list( method = method, S = S, n = n, p = p, colnames_data = colnames( data ), not.cont = not.cont, Z = Z, data = data, gcgm_NA = gcgm_NA )
    
    return( list_out )
} 
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
