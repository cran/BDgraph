## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2019  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of BDgraph package.                                                        |
#                                                                                                  |
#     BDgraph is free software: you can redistribute it and/or modify it under                     |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     This function reports below measures to assess the performance of estimated graphs:          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

compare = function( target, est, est2 = NULL, est3 = NULL, main = NULL, vis = FALSE ) 
{
    if( is.matrix( target ) ) 
    {
        if( ( sum( target == 0 ) + sum( target == 1 ) ) != ( nrow( target ) ^ 2 ) ) stop( "Element of 'target' must be 0 or 1" )
        G = target
    }
        
    if( class( target ) == "sim"   ) G <- unclass( target $ G ) 
    if( class( target ) == "graph" ) G <- unclass( target ) 
    if( ( class( target ) == "bdgraph" ) | ( class( target ) == "ssgraph" ) ) G <- BDgraph::select( target ) 
    
    if( is.matrix( est ) ) 
        if( ( sum( est == 0 ) + sum( est == 1 ) ) != ( nrow( est ) ^ 2 ) ) stop( "Element of 'est' must be 0 or 1" )
    
    if( ( class( est ) == "bdgraph" ) | ( class( est ) == "ssgraph" ) ) est <- BDgraph::select( est ) 
    if( class( est )  == "select" ) est <- est $ refit

	p = nrow( G )
	
	if( sum( dim( G ) == dim( est ) ) != 2 ) stop( " 'target' and 'est' have non-conforming size" )
	
	G[   lower.tri( G,   diag = TRUE ) ] = 0
	est[ lower.tri( est, diag = TRUE ) ] = 0
	   	   
	result        = matrix( 1, 8, 2 )
	result[ , 2 ] = compute_measures( G = G, est_G = est )
   
    if( !is.null( est2 ) )
    {
        if( is.matrix( est2 ) ) 
        {
            if( ( sum( est2 == 0 ) + sum( est2 == 1 ) ) != ( nrow( est2 ) ^ 2 ) ) stop( "Element of 'est2' must be 0 or 1" )
        }else{
            if( ( class( est2 ) == "bdgraph" ) | ( class( est2 ) == "ssgraph" ) ) est2 <- BDgraph::select( est2 ) 
            if( class( est2 ) == "select"  )  est2 <- est2 $ refit
            est2 = as.matrix( est2 )
        }
		
		if( sum( dim( G ) == dim( est2 ) ) != 2 ) stop( " 'target' and 'est2' have non-conforming size" )

		est2[ lower.tri( est2, diag = TRUE ) ] = 0

		result = cbind( result, compute_measures( G = G, est_G = est2 ) )
	} 
	
    if( !is.null( est3 ) )
    { 
        if( is.matrix( est3 ) ) 
        {
            if( ( sum( est3 == 0 ) + sum( est3 == 1 ) ) != ( nrow( est3 ) ^ 2 ) ) stop( "Element of 'est3' must be 0 or 1" )
            est3 = est3
        }else{
            if( ( class( est3 ) == "bdgraph" ) | ( class( est3 ) == "ssgraph" ) ) est3 <- BDgraph::select( est3 )
            if( class( est3 ) == "select"  )  est3 <- est3 $ refit
            est3 = as.matrix( est3 ) 
        }
		
		if( sum( dim( G ) == dim( est3 ) ) != 2 ) stop( " 'target' and 'est3' have non-conforming size" )
		
		est3[ lower.tri( est3, diag = TRUE ) ] = 0

		result = cbind( result, compute_measures( G = G, est_G = est3 ) )
	} 
	
	result[ c( 3, 4 ), 1 ]    = 0
	result[ 1        , 1 ]    = sum( G )
	result[ 2        , 1 ]    = p * ( p - 1 ) / 2 - result[ 1, 1 ]  
	result[ is.na( result ) ] = 0

	if( is.null( main ) ) 
	{
	    if( ( class( target ) == "bdgraph" ) | ( class( target ) == "ssgraph" ) )
	    {
    	    main = c( "estimate1", "estimate2" )
    		if( !is.null( est2 ) ) main = c( main, "estimate3" )
    		if( !is.null( est3 ) ) main = c( main, "estimate4" )
	    }else{
	        main = c( "True", "estimate" )
	        if( !is.null( est2 ) ) main = c( main, "estimate2" )
	        if( !is.null( est3 ) ) main = c( main, "estimate3" )
	    }
	}
		
	colnames( result ) <- main

	rownames( result ) <- c( "true positive", "true negative", "false positive", "false negative", 
                             "F1-score", "specificity", "sensitivity", "MCC" )
				
   if( vis == TRUE )
   {
		G_igraph   <- igraph::graph.adjacency( G,   mode = "undirected", diag = FALSE )
		est_igraph <- igraph::graph.adjacency( est, mode = "undirected", diag = FALSE )
		if ( p < 20 ) sizev = 15 else sizev = 2

		row_plot = ifelse( is.null( est2 ), 1, 2 )
		op       = graphics::par( mfrow = c( row_plot, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )

		igraph::plot.igraph( G_igraph,   layout = igraph::layout.circle, main = main[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		igraph::plot.igraph( est_igraph, layout = igraph::layout.circle, main = main[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		 
		if( !is.null( est2 ) )
		{
			est2_igraph <- igraph::graph.adjacency( as.matrix( est2 ), mode = "undirected", diag = FALSE )
			igraph::plot.igraph( est2_igraph, layout = igraph::layout.circle, main = main[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		if( !is.null( est3 ) )
		{ 
			est3_igraph <- igraph::graph.adjacency( as.matrix( est3 ), mode = "undirected", diag = FALSE ) 
			igraph::plot.igraph( est3_igraph, layout = igraph::layout.circle, main = main[4], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		graphics::par( op )
   }

	return( round( result, 3 ) )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    To compare measures the performance of estimated graphs based on true graph
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
compute_measures = function( G, est_G ) 
{
	upper_G     = G[     upper.tri( G     ) ]
	upper_est_G = est_G[ upper.tri( est_G ) ]
		
	tp = sum( ( upper_G == 1 ) * ( upper_est_G == 1 ) ) 
	tn = sum( ( upper_G == 0 ) * ( upper_est_G == 0 ) )
	fp = sum( ( upper_G == 0 ) * ( upper_est_G == 1 ) ) 
	fn = sum( ( upper_G == 1 ) * ( upper_est_G == 0 ) )
		
	# harmonic mean of precision and recall, called F-measure or balanced F-score
	F1score = ( 2 * tp ) / ( 2 * tp + fp + fn )

	specificity  = tn / ( tn + fp )
	sensitivity  = tp / ( tp + fn )
	# Matthews Correlation Coefficients (MCC)
	mcc          = ( ( tp * tn ) - ( fp * fn ) ) / ( sqrt( ( tp + fp ) * ( tp + fn ) ) * sqrt( ( tn + fp ) * ( tn + fn ) ) )
	
	return( c( tp, tn, fp, fn, F1score, specificity, sensitivity, mcc ) )
}
   
