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
#     To select the graph in which the edge posterior probabilities are more 
#     than "cut" value OR if cut is NULL to select the best graph ( graph with 
#     the highest posterior probability )  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

select = function( bdgraph.obj, cut = NULL, vis = FALSE )
{
	if( is.matrix( bdgraph.obj ) ) 
	{
	    if( any( bdgraph.obj < 0 ) || any( bdgraph.obj > 1 ) ) 
	        stop( "Values of matrix 'bdgraph.obj' must be between 0 and 1" )
	    
	    p_links = unclass( bdgraph.obj )
		p       = ncol( p_links )
		
	}else{
	    
	    if( ( !inherits( bdgraph.obj, "bdgraph" ) ) && ( !inherits( bdgraph.obj, "ssgraph" ) ) )
	        stop( "'bdgraph.obj' must be a matrix or an object from functions 'bdgraph()', 'bdgraph.mpl()', or 'ssgraph()'" )
	    
	    if( ( inherits( bdgraph.obj, "bdgraph" ) ) | ( inherits( bdgraph.obj, "ssgraph" ) ) )
	        p_links = bdgraph.obj $ p_links
	    
	    if( inherits( bdgraph.obj, "bdgraph" ) ) p = ncol( bdgraph.obj $ last_graph )
	    if( inherits( bdgraph.obj, "ssgraph" ) ) p = ncol( bdgraph.obj $ K_hat      )
	}
  
    if( ( is.null( p_links ) ) && ( is.null( cut ) ) )
    {
        sample_graphs <- bdgraph.obj $ sample_graphs
        graph_weights <- bdgraph.obj $ graph_weights
        
        indG_max <- sample_graphs[ which( graph_weights == max( graph_weights ) )[1] ]
        
        vec_G    <- c( rep( 0, p * ( p - 1 ) / 2 ) )
        vec_G[ which( unlist( strsplit( as.character( indG_max ), "" ) ) == 1 ) ] <- 1
        
        dimlab     <- colnames( bdgraph.obj $ last_graph )
        selected_g <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	
        selected_g[ upper.tri( selected_g ) ] <- vec_G
    
    }else{
        
        if( is.null( cut ) ) 
            cut = 0.5
        
        if( ( cut < 0 ) || ( cut > 1 ) ) 
            stop( "'cut' must be between 0 and 1" )
        
        if( is.null( p_links ) ) 
            p_links = BDgraph::plinks( bdgraph.obj, round = 10 )
        
        selected_g                   = 0 * p_links
        selected_g[ p_links >  cut ] = 1
        selected_g[ p_links <= cut ] = 0
    }
    
	if( vis )
	{
	    main = ifelse( is.null( cut ), "Graph with highest posterior probability.", 
	                   paste( c( "Graph with links posterior probabilities > ",  cut ), collapse = "" ) )
	    
	    BDgraph::plot.graph( selected_g, main = main )
	}

	return( selected_g )
}
       
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
