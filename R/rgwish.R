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
#     Sampling from G-Wishart distribution                                     |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

rgwish = function( n = 1, adj = NULL, b = 3, D = NULL, threshold = 1e-8 )
{
	if( b <= 2 )         stop( "'b' must be more than 2" )
	if( is.null( adj ) ) stop( "'adj' must be determined" )

    if( is.matrix( adj )           ) G <- unclass( adj )
  # if( inherits( adj, "graph" )   ) G <- unclass( adj )
    if( inherits( adj, "sim" )     ) G <- unclass( adj $ G )
    if( inherits( adj, "bdgraph" ) ) G <- BDgraph::select( adj ) 
    if( inherits( adj, "ssgraph" ) ) G <- BDgraph::select( adj ) 
    
    
    if( ( sum( G == 0 ) + sum( G == 1 ) ) != ( nrow( G ) ^ 2 ) ) 
		stop( "Elements of matrix 'adj' must be 0 or 1" )
    
    G <- as.matrix( G )
    diag( G ) <- 0
    
    if( !isSymmetric( G ) )
    {
        G[ lower.tri( G ) ] <- 0
        G                   <- G + t( G )
    }
	
	p <- nrow( G )
	if( p < 1 ) stop( "'p' must be more than or equal with 1" )
	
	if( is.null( D )      ) D <- diag( p )
	if( !isSymmetric( D ) ) stop( "'D' must be a positive definite matrix" )
	if( nrow( D ) != p    ) stop( "'G' and 'D' dimentions differ" )
	
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
			result = .C( "rgwish_c", as.integer(G), as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
			samples[ , , i ] = matrix( result $ K, p, p ) 		
		}
	}else{
	
		result  = .C( "rgwish_c", as.integer(G), as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
		samples = matrix( result $ K, p, p ) 		
	}

	return( samples )   
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

