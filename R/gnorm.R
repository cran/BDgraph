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
#     To compute the normalizing constant of G-Wishart distribution based on   |
#     Monte Carlo algorithm according to below paper                           |   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Atay-Kayis & Massam (2005). A monte carlo method for computing the       |
#  marginal likelihood in nondecomposable Gaussian graphical models, Biometrika|    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

gnorm = function( adj, b = 3, D = diag( ncol( adj ) ), iter = 100 )
{
	if ( b < 3 )         stop( "'b' must be more than 2" )
	if( is.null( adj ) ) stop( "'adj' must be determined" )

    G <- unclass( adj )
    G <- as.matrix( G )
    p <- nrow( G )

    if( p != ncol( G ) ) stop( "'adj' must be a square matrix" )
    if( ( sum( G == 0 ) + sum( G == 1 ) ) != ( nrow( G ) ^ 2 ) ) stop( "Elements of matrix 'adj' must be 0 or 1" )

	G[ lower.tri( G, diag = TRUE ) ] <- 0

	Ti      = chol( solve( D ) )
	H       = Ti / t( matrix( rep( diag( Ti ), p ), p, p ) )
	check_H = identical( H, diag( p ) ) * 1

	nu         = rowSums( G )
	size_graph = sum( G )
    	
	# For the case, G is a full graph
	if( size_graph == ( p * ( p - 1 ) / 2 ) )
	{
	    logIg = ( size_graph / 2 ) * log( pi ) + ( p * ( b + p - 1 ) / 2 ) * log( 2 ) +
	        sum( lgamma( ( b + nu ) / 2 ) ) - ( ( b + p - 1 ) / 2 ) * log( det( D ) )
	}
	
	# For the case, G is an empty graph
	if( size_graph == 0 )
	    logIg = ( p * b / 2 ) * log( 2 ) + p * lgamma( b / 2 ) - ( b / 2 ) * sum( log( diag( D ) ) )
	    
	if( ( size_graph != ( p * ( p - 1 ) / 2 ) ) & ( size_graph != 0 ) )
	{   # For the case G is NOT full graph
    	# - -  Monte Carlo glorithm which is implemented in C++ - - - - - - - -|
    	f_T = c( rep( 0, iter ) )
    	result = .C( "log_exp_mc", as.integer( G ), as.integer( nu ), as.integer( b ), as.double( H ), as.integer( check_H ), 
    	              as.integer( iter ), as.integer( p ), f_T = as.double( f_T ), PACKAGE = "BDgraph" )
    	
    	f_T      = c( result $ f_T )
    	log_Ef_T = log( mean( exp( - f_T / 2 ) ) )
    	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    	
    	c_dT = ( size_graph / 2 ) * log( pi ) + ( p * b / 2 + size_graph ) * log( 2 ) +
    	    sum( lgamma( ( b + nu ) / 2 ) ) + sum( ( b + nu + colSums( G ) ) * log( diag( Ti ) ) )
    	
    	logIg = c_dT + log_Ef_T
	}
	
	return( logIg ) 
}
  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

