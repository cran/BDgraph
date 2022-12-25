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

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#  Get graph from class objects "sim", "graph", "bdgraph", "ssgraph", or "select"
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_graph = function( obj_G, cut = 0.5 )
{
    if( is.matrix( obj_G ) ) 
    {
        if( nrow( obj_G ) != ncol( obj_G ) ) stop( "Adjacency matrix must be squere" )
        
        if( ( sum( obj_G == 0 ) + sum( obj_G == 1 ) ) != ( ncol( obj_G ) ^ 2 ) ) 
            stop( "Elements of adjacency matrix must be 0 or 1" )
        
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
        if( ( g.prior <= 0 ) | ( g.prior >= 1 ) ) stop( "'g.prior' must be between 0 and 1" )
        g.prior = matrix( g.prior, p, p )
    }else{
        if( ( nrow( g.prior ) != p ) | ( ncol( g.prior ) != p ) ) stop( "'g.prior' and 'data' have non-conforming size" )
        if( any( g.prior < 0 ) || any( g.prior > 1 ) ) stop( "Elements of  matrix 'g.prior' must be between 0 and 1" )
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
        if( ( sum( g.start == 0 ) + sum( g.start == 1 ) ) != ( p ^ 2 ) ) 
            stop( "Elements of matrix 'g.start' must be 0 or 1" )
        
        G = g.start
    }
    
    if( inherits( g.start, "sim"   ) ) G <- unclass( g.start $ G )
    if( inherits( g.start, "graph" ) ) G <- unclass( g.start )
    
    if( ( inherits( g.start, "bdgraph"   ) ) |  ( inherits( g.start, "ssgraph" ) ) ) G <- g.start $ last_graph
    if( ( inherits( g.start, "character" ) ) && ( g.start == "empty" ) ) G = matrix( 0, p, p )
    if( ( inherits( g.start, "character" ) ) && ( g.start == "full"  ) ) G = matrix( 1, p, p )
    
    if( ( nrow( G ) != p ) | ( ncol( G ) != p ) ) 
        stop( "'g.start' and 'data' have non-conforming size" )
    
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
        K      = matrix( result $ K, p, p ) 
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
    
    if( !is.matrix( data ) & !is.data.frame( data ) ) 
        stop( "'data' must be a matrix or dataframe" )
    
    if( is.data.frame( data ) ) 
        data <- data.matrix( data )
    
    if( any( is.na( data ) ) ) 
    {
        if( method == "dw"  ) stop( "'bdgraph.dw()' does not deal with missing values" )	
        if( method == "ggm" ) stop( "'ggm' method does not deal with missing values. You could choose option 'method = \"gcgm\"'" )	
        if( method == "tgm" ) stop( "'tgm' method does not deal with missing values. You could choose option 'method = \"gcgm\"'" )	
        
        gcgm_NA = 1
    }else{
        gcgm_NA = 0
    }
    
    if( isSymmetric( data ) )
    {
        if( method == "gcgm" ) stop( "'method = \"gcgm\"' requires all data" )
        if( method == "dw"   ) stop( "'method = \"dw\"' requires all data" )
        if( method == "tgm"  ) stop( "'method = \"tgm\"' requires all data" )
        if( is.null( n )     ) stop( "Please specify the number of observations 'n'" )
    }
    
    p <- ncol( data )
    
    if( p < 3 ) 
        stop( "Number of variables/nodes ('p') must be more than 2" )
    
    if( is.null( n ) ) 
        n <- nrow( data )
    
    if( method == "ggm" ) 
    {
        if( isSymmetric( data ) )
        {
            cat( "Input is identified as the covariance matrix \n" )
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
                if( length( unique( data[ , j ] ) ) > min( n / 2 ) ) 
                    not.cont[ j ] = 0
        }else{
            
            if( !is.vector( not.cont )  ) 
                stop( "'not.cont' must be a vector with length of number of variables" )
            
            if( length( not.cont ) != p ) 
                stop( "'not.cont' must be a vector with length of number of variables" )
            
            if( ( sum( not.cont == 0 ) + sum( not.cont == 1 ) ) != p ) 
                stop( "Elements of vector 'not.cont' must be 0 or 1" )
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
            if( !is.vector( not.cont )  ) stop( "'not.cont' must be a vector with length of number of variables" )
            if( length( not.cont ) != p ) stop( "'not.cont' must be a vector with length of number of variables" )
            if( ( sum( not.cont == 0 ) + sum( not.cont == 1 ) ) != p ) stop( "Elements of vector 'not.cont' must be 0 or 1" )
        }
        
        # for non-Gaussian data
        Z                  <- stats::qnorm( apply( data, 2, rank, ties.method = "random" ) / ( n + 1 ) )
        Zfill              <- matrix( stats::rnorm( n * p ), n, p )   # for missing values
        Z[ is.na( data ) ] <- Zfill[ is.na( data ) ]                  # for missing values
        Z                  <- t( ( t( Z ) - apply( Z, 2, mean ) ) / apply( Z, 2, stats::sd ) )
        S                  <- t( Z ) %*% Z
    } 
    
    if( method == "tgm" )
        S <- t( data ) %*% data
    
    list_out = list( method = method, S = S, n = n, p = p, colnames_data = colnames( data ), data = data )
    
    if( method == "gcgm" )
    {
        list_out$not.cont = not.cont
        list_out$gcgm_NA  = gcgm_NA
        list_out$Z        = Z
        list_out$R        = R
    }

    if( method == "dw" )
    {
        list_out$not.cont = not.cont
        list_out$gcgm_NA  = gcgm_NA
        list_out$Z        = Z
    }
        
    return( list_out )
} 
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_Ds_tgm_R = function( data, tu, mu, D, n, p, in_C = TRUE )
{
    if( in_C == TRUE )
    {
        S  = matrix( 0, p, p )
        Ds = matrix( 0, p, p )
        
        # void get_Ds_tgm( double data[], double D[], double mu[], double tu[], double Ds[], double S[], int *n, int *p )
        result = .C( "get_Ds_tgm", as.double(data), as.double(D), as.double (mu), as.double(tu), 
                     Ds = as.double(Ds), S = as.double(S), as.integer(n), as.integer(p), PACKAGE = "BDgraph" )
        
        S  = matrix( result $ S , p, p ) 
        Ds = matrix( result $ Ds, p, p ) 
    }else{
        
        #X_new = matrix( nrow = n, ncol = p )
    	#for( i in 1:n ){
        #    sqrt_tu_i = sqrt( tu[ i ] );
        #    for( j in 1:p ) X_new[ i, j ] = sqrt_tu_i * ( data[ i, j ] - mu[ j ] );
    	#}
    	#S  = t( X_new ) %*% X_new
    	
    	S  = matrix( 0, p, p )
    	for( i in 1:p )
    	    for( j in 1:p )
    	        for( k in 1:n )
    	            S[ i, j ] = S[ i, j ] + tu[ k ] * ( data[ k, i ] - mu[ i ] ) * ( data[ k, j ] - mu[ j ] );
        
    	Ds = D + S
    }
    
	return( list( S = S, Ds = Ds ) )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
get_Ts_R = function( Ds, in_C = TRUE )
{
    if( in_C == TRUE )
    {
        p = ncol( Ds )
        
        Ts      = matrix( 0, p, p )
        inv_Ds  = matrix( 0, p, p )
        copy_Ds = matrix( 0, p, p )
       
        # void get_Ts( double Ds[], double Ts[], double inv_Ds[], double copy_Ds[], int *p )
        result = .C( "get_Ts", as.double(Ds), Ts = as.double(Ts), as.double(inv_Ds), as.double(copy_Ds), as.integer(p), PACKAGE = "BDgraph" )
        
        Ts = matrix( result $ Ts, p, p ) 
        
    }else{
    	invDs = solve( Ds )
	    Ts    = chol( invDs )
    }

	return( Ts )    
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
update_tu_R = function( tu, data, K, mu, nu, n, p, in_C = TRUE )
{
    if( in_C == TRUE )
    {
        Ts      = matrix( 0, p, p )
        inv_Ds  = matrix( 0, p, p )
        copy_Ds = matrix( 0, p, p )
       
        # void update_tu( double data[], double K[], double tu[], double mu[], double *nu, int *n, int *p )
        result = .C( "update_tu", as.double(data), as.double(K), tu = as.double(tu), as.double(mu), as.double(nu), as.integer(n), as.integer(p), PACKAGE = "BDgraph" )
        
        tu = c( result $ tu ) 
        
    }else{
        d_mu_i = numeric( p )
        
        for( i in 1:n )
    	{
            for( j in 1:p )
                d_mu_i[ j ] = data[ i, j ] - mu[ j ];
                
            #d_mu_i_x_K = d_mu_i %*% K;
            #delta_y_i  = sum( d_mu_i_x_K * d_mu_i );
            
            delta_y_i = 0.0;
            for( k in 1:p )
                for( l in 1:p )
                    delta_y_i = delta_y_i + d_mu_i[ l ] * K[ l, k ] * d_mu_i[ k ];
                
    		shape_tu_i = ( nu + p ) / 2.0;
    		rate_tu_i  = ( nu + delta_y_i ) / 2.0;
    			
    		tu[ i ] = stats::rgamma( 1, shape = shape_tu_i, scale = 1.0 / rate_tu_i )
        }
    }

	return( tu )    
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
update_mu_R = function( data, mu, tu, n, p, in_C = TRUE )
{
    if( in_C == TRUE )
    {
        # void update_mu( double data[], double mu[], double tu[], int *n, int *p )
        result = .C( "update_mu", as.double(data), mu = as.double(mu), as.double(tu), as.integer(n), as.integer(p), PACKAGE = "BDgraph" )
        
        mu = c( result $ mu ) 
        
    }else{
        mu = c( tu %*% data / sum( tu ) )
    }

	return( mu )    
}
                      
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |



















