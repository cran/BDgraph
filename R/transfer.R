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
#     For count data: to transfer raw data for the algorithm                   |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

transfer = function( r_data )
{
	if( inherits( r_data, "sim" ) ) r_data <- r_data $ data
  
	n    = dim( r_data )[ 1 ]
	p    = dim( r_data )[ 2 ]
	data = matrix( 0, nrow = n, ncol = p + 1 )
	size_unique_data = 0
	
	result = .C( "transfer_data", 
				 as.integer(r_data), 
				 data = as.integer(data), 
				 as.integer(n), 
				 as.integer(p), 
				 size_unique_data = as.integer(size_unique_data), 
				 PACKAGE = "BDgraph" )
	
	size_unique_data = result $ size_unique_data
	
	label = colnames( r_data )
	if( is.null( label ) ) label = 1:p
  
	data = matrix ( result $ data, n, p + 1, dimnames = list( NULL, c( label, "ferq" ) ) ) 
	data = data[ (1:size_unique_data), ]

	return( data )
}
  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

