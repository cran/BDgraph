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
#     To compute the normalizing constant of G-Wishart distribution based on Monte Carlo           |
#     algorithm according to below paper                                                           |
## ------------------------------------------------------------------------------------------------|
#     Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the marginal         |
#     likelihood in nondecomposable Gaussian graphical models, Biometrika, 92(2):317-335           |
## ------------------------------------------------------------------------------------------------|

gnorm = function( adj.g, b = 3, D = diag( ncol( adj.g ) ), iter = 100 )
{
	if ( b < 3 )           stop( "Parameter 'b' must be more than 2" )
	if( is.null( adj.g ) ) stop( "Adjacency matrix should be determined" )

    G <- unclass( adj.g )
    G <- as.matrix( G )
	if( sum( ( G == 1 ) * ( G == 0 ) ) != 0 ) stop( "Elements of matrix G must be zero or one" )	

	G[ lower.tri( G, diag = TRUE ) ] <- 0

	p       = nrow( G )
	Ti      = chol( solve( D ) )
	H       = Ti / t( matrix( rep( diag( Ti ), p ), p, p ) )
	check_H = identical( H, diag( p ) ) * 1

	nu  = rowSums( G )
	f_T = c( rep( 0, iter ) )

## ---- main BDMCMC algorithms implemented in C++ -------------------------------------------------|
	result = .C( "log_exp_mc", as.integer( G ), as.integer( nu ), as.integer( b ), as.double( H ), as.integer( check_H ), 
	              as.integer( iter ), as.integer( p ), f_T = as.double( f_T ), PACKAGE = "BDgraph" )
	f_T    = c( result $ f_T )
## ------------------------------------------------------------------------------------------------|
	
	log_Ef_T = log( mean( exp( - f_T / 2 ) ) )

	size_graph = sum( G )
	c_dT = ( size_graph / 2 ) * log( pi ) + ( p * b / 2 + size_graph ) * log( 2 ) +
	          sum( lgamma( ( b + nu ) / 2 ) ) + sum( ( b + nu + colSums( G ) ) * log( diag( Ti ) ) )
	  
	logIg = c_dT + log_Ef_T
	
	return( logIg ) 
}
      
