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
#     Sampling from Wishart distribution                                       |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

rwish = function( n = 1, p = 2, b = 3, D = diag( p ) )
{
    if( p < 1             ) stop( "'p' must be more than or equal with 1" )
    if( b <= 2            ) stop( "'b' must be more than 2" )
	if( !isSymmetric( D ) ) stop( "'D' must be a positive definite matrix" )
    if( n < 1             ) stop( "'n' must be more than or equal with 1" )
    if( ncol( D ) != p    ) stop( "'p' and 'D' have non-conforming size" )
    
	Ti = chol( solve( D ) ) 
	K  = matrix( 0, p, p )

	if( n > 1 )
	{
		samples = array( 0, c( p, p, n ) )
		
		for ( i in 1 : n )
		{
			result       = .C( "rwish_c", as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
			samples[,,i] = matrix( result $ K, p, p ) 		
		}
	}else{
		result  = .C( "rwish_c", as.double(Ti), K = as.double(K), as.integer(b), as.integer(p), PACKAGE = "BDgraph" )
		samples = matrix( result $ K, p, p ) 				
	}	

	return( samples )   
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

