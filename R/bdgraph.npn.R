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
#     Non-parametric transfer function for non-normal data                     |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bdgraph.npn = function( data, npn = "shrinkage", npn.thresh = NULL )
{
    if( inherits( data, "sim" ) ) data <- data $ data
    
    if( !is.matrix( data ) & !is.data.frame( data ) ) stop( "Data must be a matrix or dataframe" )	
    if( is.data.frame( data ) ) data = data.matrix( data )
    if( any( is.na( data ) ) ) stop( "Data must contain no missing values" ) 
	
	n <- nrow( data )
	
  	# - - - shrinkage transfer - - - - - - - - - - - - - - - - - - - - - - - - |
	
	if( npn == "shrinkage" )
	{
		data = stats::qnorm( apply( data, 2, rank ) / ( n + 1 ) )
		#data = data / stats::sd( data[ , 1 ] )
		data = t( ( t( data ) - apply( data, 2, mean ) ) / apply( data, 2, stats::sd ) )
	}
	
	# - - - truncation transfer - - - - - - - - - - - - - - - - - - - - - - - -|
	
	if( npn == "truncation" )
	{
		if( is.null( npn.thresh ) ) npn.thresh = 0.25 * ( n ^ ( -0.25 ) ) * ( pi * log( n ) ) ^ ( -0.5 )
		
		data = stats::qnorm( pmin( pmax( apply( data, 2, rank ) / n, npn.thresh ), 1 - npn.thresh ) )
    	data = data / stats::sd( data[ , 1 ] )
	}

	# - - - skeptic transfer - - - - - - - - - - - - - - - - - - - - - - - - - |
	
	if( npn == "skeptic" ) data = 2 * sin( pi / 6 * stats::cor( data, method = "spearman" ) )
	
	return( data )
}
 
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

