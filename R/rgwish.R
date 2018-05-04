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
#     Sampling from G-Wishart distribution                                                         |
## ------------------------------------------------------------------------------------------------|

rgwish = function( n = 1, adj.g = NULL, b = 3, D = NULL )
{
	if( b <= 2 )           stop( "For G-Wishart distribution parameter 'b' must be more than 2" )
	if( is.null( adj.g ) ) stop( "Adjacency matrix must be determined" )

	if( class( adj.g ) == "sim"   ) G <- adj.g $ G
	if( class( adj.g ) == "graph" ) G <- unclass( adj.g )
	
    G <- as.matrix( G )
    
	if( sum( ( G == 1 ) * ( G == 0 ) ) != 0 ) stop( "Elements of matrix G must be 0 or 1" )	

	if( !isSymmetric( G ) )
	{
		G[ lower.tri( G, diag( TRUE ) ) ] <- 0
		G  = G + t( G )
	}
	
	diag( G ) = 0
	
	p <- nrow( G )
	if( p < 1 ) stop( "'p' must be more than or equal with 1" )
	
	if( is.null( D ) ) D <- diag( p )
	if( !isSymmetric( D ) ) stop( "Matrix 'D' must be positive definite matrix." )
		
	if( dim( D )[ 1 ] != p ) stop( "Dimension of matrix G and D must be the same." )
	
	if( p == 1 )
	    return( rwish( n = n, p = p, b = b, D = D ) )

	if( sum( G ) == ( p * ( p - 1 ) ) )
	    return( rwish( n = n, p = p, b = b, D = D ) )
	    
	Ti = chol( solve( D ) )
	K  = matrix( 0, p, p )
	
	if( n > 1 )
	{
		samples = array( 0, c( p, p, n ) )
		
		for( i in 1 : n )
		{
			result          = .C( "rgwish_c", as.integer(G), as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
			samples[ , , i] = matrix( result $ K, p, p ) 		
		}
	}else{
	
		result  = .C( "rgwish_c", as.integer(G), as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
		samples = matrix( result $ K, p, p ) 		
	}

	return( samples )   
}
  
