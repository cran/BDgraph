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
#     Sampling from COMPLEX G-Wishart distribution                                                 |
## ------------------------------------------------------------------------------------------------|

rgcwish = function( n = 1, adj.g = NULL, b = 3, D = NULL )
{
	if( b <= 2 )           stop( "For complex G-Wishart distribution parameter 'b' must be more than 2" )
	if( is.null( adj.g ) ) stop( "Adjacency matrix must be determined" )
    
    if( is.matrix( adj.g )          ) G <- unclass( adj.g )
    #if( class( adj.g ) == "graph"   ) G <- unclass( adj.g )
    if( class( adj.g ) == "sim"     ) G <- adj.g $ G
    if( class( adj.g ) == "bdgraph" ) G <- BDgraph::select( adj.g ) 
    if( class( adj.g ) == "ssgraph" ) G <- BDgraph::select( adj.g ) 
    
    if( !isSymmetric( G ) )
    {
        G[ lower.tri( G ) ] <- 0
        G  = G + t( G )
    }
    
    
    if( sum( ( G == 1 ) * ( G == 0 ) ) != 0 ) stop( "Elements of matrix 'adj.g' must be 0 or 1" )	

	G <- as.matrix( G )
	diag( G ) = 0
	
	p <- nrow( G )  
	if( p < 1 ) stop( "'p' must be more than or equal with 1" )
	
	if( is.null( D ) ) D <- diag( p )
	if( !isSymmetric( D ) ) stop( "Matrix 'D' must be positive definite" )
	
	if( nrow( D ) != p ) stop( "Dimension of matrix G and D must be the same." )
		
	inv_D = solve( D )
	row1  = cbind( Re(inv_D), -Im(inv_D) )
	row2  = cbind( Im(inv_D), Re(inv_D) )
	sig   = rbind( row1, row2 ) / 2
	Ls    = t( chol(sig) )
	K     = matrix( 0, p, p )
	
	if( n > 1 )
	{
		samples = array( 0, c( p, p, n ) )
		for( i in 1 : n )
		{
			result       = .C( "rgcwish_c", as.integer(G), as.double(Ls), K = as.complex(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
			samples[,,i] = matrix( result $ K, p, p ) 		
		}
	}else{
		result  = .C( "rgcwish_c", as.integer(G), as.double(Ls), K = as.complex(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
		samples = matrix( result $ K, p, p ) 				
	}	

	return( samples )   
}
  
