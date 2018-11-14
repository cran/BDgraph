## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2018  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of BDgraph package.                                                        |
#                                                                                                  |
#     BDgraph is free software: you can redistribute it and/or modify it under                     |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Plot of graph size to check the convergency of BDMCMC algorithm                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

traceplot = function( bdgraph.obj, acf = FALSE, pacf = FALSE, main = NULL, ... )
{
    if( ( class( bdgraph.obj ) == "bdgraph" ) | ( class( bdgraph.obj ) == "ssgraph" ) )
    {
        if( is.null( bdgraph.obj $ all_graphs ) ) stop( "'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option save = TRUE" )
        if( is.null( bdgraph.obj $ all_graphs ) ) stop( "'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option save = TRUE" )
    }else{
        stop( "'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()'" )
    }
    
	sample_graphs     = bdgraph.obj $ sample_graphs
    all_graphs        = bdgraph.obj $ all_graphs
	graph_weights     = bdgraph.obj $ graph_weights
	sizesample_graphs = sapply( sample_graphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )  
	sizeall_graphs    = sizesample_graphs[ all_graphs ]
	which_G_max       = which( max( graph_weights ) == graph_weights )
	size_selected_g   = sizeall_graphs[ which_G_max ] 

	if( is.null( main ) ) main = "Trace of graph size"
	
	x_vec = 1 : length( all_graphs )
	
	if ( acf == FALSE & pacf == FALSE )
	{
		graphics::plot( x = x_vec, sizeall_graphs, type = "l", main = main, cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, ylab = "Graph size", xlab = "Iteration", ... )
		graphics::abline( h = size_selected_g, col = "red" )	   
	}
	
	if ( acf == TRUE & pacf == TRUE )
	{
		op = graphics::par( mfrow = c( 2, 2 ), pty = "s" )  
		graphics::plot( x = x_vec, sizeall_graphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		graphics::abline( h = size_selected_g, col = "red" )	  
		acf( sizeall_graphs,  main = "ACF for graph size" )
		pacf( sizeall_graphs, main = "PACF for graph size" )
		graphics::par( op )
	}
	
	if ( acf == TRUE & pacf == FALSE )
	{
		op <- graphics::par( mfrow = c( 1, 2 ), pty = "s" ) 
		graphics::plot( x = x_vec, sizeall_graphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		graphics::abline( h = size_selected_g, col = "red" )	  
		acf( sizeall_graphs, main = "ACF for graph size" )
		graphics::par( op )
	}
	
	if ( acf == FALSE & pacf == TRUE )
	{
		op <- graphics::par( mfrow = c( 1, 2 ), pty = "s" ) 
		graphics::plot( x = x_vec, sizeall_graphs, type = "l", main = main, ylab = "Graph size", xlab = "Iteration", ... )
		graphics::abline( h = size_selected_g, col = "red" )	  
		pacf( sizeall_graphs, main = "PAIC for graph size" )
		graphics::par( op )
	}		
}  
      
