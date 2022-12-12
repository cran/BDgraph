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
#     Computing posterior probabilities of all possible links                  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

plinks = function( bdgraph.obj, round = 2, burnin = NULL )
{
    if( ( !inherits( bdgraph.obj, "bdgraph" ) ) && ( !inherits( bdgraph.obj, "ssgraph" ) ) )
        stop( "'bdgraph.obj' must be an object of functions 'bdgraph()', 'bdgraph.mpl()', 'bdgraph.dw()', or 'ssgraph()'" )
    
    if( inherits( bdgraph.obj, "bdgraph" ) )
    {
        if( is.null( bdgraph.obj $ sample_graphs ) )
    	  {
    		    p_links = bdgraph.obj $ p_links
    	  }else{
    		
        		p     <- nrow( bdgraph.obj $ last_graph )
        		vec_G <- numeric( length = ( p * ( p - 1 ) / 2 ) )
        
        		sample_graphs <- bdgraph.obj $ sample_graphs
        		graph_weights <- bdgraph.obj $ graph_weights
        
        		if( is.null( burnin ) )
        		{
          			for( i in 1 : length( sample_graphs ) )
          			{
            				inp          <- which( unlist( strsplit( as.character( sample_graphs[ i ] ), "" ) ) == 1 )
            				vec_G[ inp ] <- vec_G[ inp ] + graph_weights[ i ]
          			}
          			
          			sum_graph_weights = sum( graph_weights )
          	}else{
    
          			all_graphs  <- bdgraph.obj $ all_graphs
          			all_weights <- bdgraph.obj $ all_weights
          			
          			sum_graph_weights <- 0
          		   
          			for( i in ( burnin + 1 ) : length( all_graphs ) )
          			{
            				inp               <- which( unlist( strsplit( as.character( sample_graphs[ all_graphs[ i ] ] ), "" ) ) == 1 )
            				vec_G[ inp ]      <- vec_G[ inp ] + all_weights[ i ]
            				sum_graph_weights <- sum_graph_weights + all_weights[ i ]
          			}
    		    }
    		
        		label   <- colnames( bdgraph.obj $ last_graph ) 
        		p_links <- matrix( 0, p, p, dimnames = list( label, label ) )
        		p_links[ upper.tri( p_links ) ] <- vec_G / sum_graph_weights	
    	  }
    }

    if( inherits( bdgraph.obj, "ssgraph" ) )
        p_links = bdgraph.obj $ p_links        
        
	  return( round( p_links, round ) )
}
          
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
