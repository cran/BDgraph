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
#     Plot of graph size to check the convergency of BDMCMC algorithm          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

traceplot = function( bdgraph.obj, acf = FALSE, pacf = FALSE, main = NULL, ... )
{
    if( ( inherits( bdgraph.obj, "bdgraph" ) ) | ( inherits( bdgraph.obj, "ssgraph" ) ) )
    {
        if( is.null( bdgraph.obj $ all_graphs ) ) 
            stop( "'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option 'save = TRUE'" )

    	sample_graphs     = bdgraph.obj $ sample_graphs
        all_graphs        = bdgraph.obj $ all_graphs
    	graph_weights     = bdgraph.obj $ graph_weights
    	
    	sizesample_graphs = sapply( sample_graphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )  
    	
    	sizeall_graphs    = sizesample_graphs[ all_graphs ]
    	which_G_max       = which( max( graph_weights ) == graph_weights )
    	size_selected_g   = sizeall_graphs[ which_G_max ] 
    	
    	sample_mcmc = sizeall_graphs
    
    	if( is.null( main ) ) main = "Trace of graph size"
    	ylab = "Graph size"
    
    }else{
        
        if( !is.vector( bdgraph.obj ) )
            stop( "'bdgraph.obj' must be an object of functions 'bdgraph()', 'bdgraph.mpl()', or 'ssgraph()' or a vector" )
        
        sample_mcmc = bdgraph.obj
        
        if( is.null( main ) ) main = "Trace of MCMC sample"
        ylab = ""
    }
    
	if( acf == FALSE & pacf == FALSE ) op = graphics::par( mfrow = c( 1, 1 ), pty = "s" )
	if( acf == TRUE  & pacf == TRUE  ) op = graphics::par( mfrow = c( 2, 2 ), pty = "s" ) 
	if( acf == TRUE  & pacf == FALSE ) op = graphics::par( mfrow = c( 1, 2 ), pty = "s" )
	if( acf == FALSE & pacf == TRUE  ) op = graphics::par( mfrow = c( 1, 2 ), pty = "s" )
    
	x_vec = 1 : length( sample_mcmc )
	
	graphics::plot( x = x_vec, 
	                y = sample_mcmc, 
	                type = "l", 
	                main = main, 
	                col = "lightblue", 
	                ylab = ylab, 
	                xlab = "Iteration", ... )

	if( !is.vector( bdgraph.obj ) )	
	    graphics::lines( x = x_vec, y = rep( size_selected_g, length( sample_mcmc ) ), 
	                     col = "blue" )
	
	if( acf  == TRUE ) 
	    acf(  sample_mcmc, main = "ACF for graph size" )
	
	if( pacf == TRUE ) 
	    pacf( sample_mcmc, main = "PACF for graph size" )
	
	graphics::par( op )
}  
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
