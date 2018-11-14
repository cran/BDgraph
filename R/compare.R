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
#     This function reports below measures to assess the performance of estimated graphs:          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# True positive:  number of correctly estimated links.                                             |
# True negative:  number of true non-existing links which is correctly estimated.                  |
# False positive: number of links which are not in the true graph, but are incorrectly estimated.  |
# False negative: number of links which they are in the true graph, but are not estimated.         |
# F1-score:       weighted average of the positive predictive and true positive rate.              |
# Specificity:    Specificity value reaches.                                                       |
# Sensitivity:    Sensitivity value reaches.                                                       |
# MCC:            Matthews Correlation Coefficients (MCC).                                         |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

compare = function( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL, 
                    main = NULL, vis = FALSE ) 
{
    if( is.matrix( sim.obj ) ) 
    {
        if( ( sum( sim.obj == 0 ) + sum( sim.obj == 1 ) ) != ( nrow( sim.obj ) ^ 2 ) ) stop( "Element of 'sim.obj' must be 0 or 1" )
        G = sim.obj
    }
        
    if( class( sim.obj ) == "sim"   ) G <- unclass( sim.obj $ G ) 
    if( class( sim.obj ) == "graph" ) G <- unclass( sim.obj ) 
    if( ( class( sim.obj ) == "bdgraph" ) | ( class( sim.obj ) == "ssgraph" ) )
    {
        G <- BDgraph::select( sim.obj ) 
    }
    
    if( is.matrix( bdgraph.obj ) ) 
    {
        if( ( sum( bdgraph.obj == 0 ) + sum( bdgraph.obj == 1 ) ) != ( nrow( bdgraph.obj ) ^ 2 ) ) stop( "Element of 'bdgraph.obj' must be 0 or 1" )
        est = bdgraph.obj
    }
    
    if( ( class( bdgraph.obj ) == "bdgraph" ) | ( class( bdgraph.obj ) == "ssgraph" ) ) est <- BDgraph::select( bdgraph.obj ) 
    if( class( bdgraph.obj )  == "select" ) est <- bdgraph.obj $ refit

	p = nrow( G )
	
	if( sum( dim( G ) == dim( est ) ) != 2 ) stop( "'sim.obj' and 'bdgraph.obj' have non-conforming size" )
	
	G[   lower.tri( G,   diag = TRUE ) ] = 0
	est[ lower.tri( est, diag = TRUE ) ] = 0
	   	   
	result        = matrix( 1, 8, 2 )
	result[ , 2 ] = compute_measures( G = G, est_G = est )
   
    if( !is.null( bdgraph.obj2 ) )
    {
        if( is.matrix( bdgraph.obj2 ) ) 
        {
            if( ( sum( bdgraph.obj2 == 0 ) + sum( bdgraph.obj2 == 1 ) ) != ( nrow( bdgraph.obj2 ) ^ 2 ) ) stop( "Element of 'bdgraph.obj2' must be 0 or 1" )
            est2 = bdgraph.obj2
        }else{
            if( ( class( bdgraph.obj2 ) == "bdgraph" ) | ( class( bdgraph.obj2 ) == "ssgraph" ) ) est2 <- BDgraph::select( bdgraph.obj2 ) 
            if( class( bdgraph.obj2 ) == "select"  )  est2 <- bdgraph.obj2 $ refit
            est2 = as.matrix( est2 )
        }
		
		if( sum( dim( G ) == dim( est2 ) ) != 2 ) stop( " 'sim.obj' and 'bdgraph.obj2' have non-conforming size" )

		est2[ lower.tri( est2, diag = TRUE ) ] = 0

		result = cbind( result, compute_measures( G = G, est_G = est2 ) )
	} 
	
    if( !is.null( bdgraph.obj3 ) )
    { 
        if( is.matrix( bdgraph.obj3 ) ) 
        {
            if( ( sum( bdgraph.obj3 == 0 ) + sum( bdgraph.obj3 == 1 ) ) != ( nrow( bdgraph.obj3 ) ^ 2 ) ) stop( "Element of 'bdgraph.obj3' must be 0 or 1" )
            est3 = bdgraph.obj3
        }else{
            if( ( class( bdgraph.obj3 ) == "bdgraph" ) | ( class( bdgraph.obj3 ) == "ssgraph" ) ) est3 <- BDgraph::select( bdgraph.obj3 )
            if( class( bdgraph.obj3 ) == "select"  )  est3 <- bdgraph.obj3 $ refit
            est3 = as.matrix( est3 ) 
        }
		
		if( sum( dim( G ) == dim( est3 ) ) != 2 ) stop( "'sim.obj' and 'bdgraph.obj3' have non-conforming size" )
		
		est3[ lower.tri( est3, diag = TRUE ) ] = 0

		result = cbind( result, compute_measures( G = G, est_G = est3 ) )
	} 
	
	result[ c( 3, 4 ), 1 ]    = 0
	result[ 1        , 1 ]    = sum( G )
	result[ 2        , 1 ]    = p * ( p - 1 ) / 2 - result[ 1, 1 ]  
	result[ is.na( result ) ] = 0

	if( is.null( main ) ) 
	{
	    if( ( class( sim.obj ) == "bdgraph" ) | ( class( sim.obj ) == "ssgraph" ) )
	    {
    	    main = c( "estimate1", "estimate2" )
    		if( !is.null( bdgraph.obj2 ) ) main = c( main, "estimate3" )
    		if( !is.null( bdgraph.obj3 ) ) main = c( main, "estimate4" )
	    }else{
	        main = c( "True", "estimate" )
	        if( !is.null( bdgraph.obj2 ) ) main = c( main, "estimate2" )
	        if( !is.null( bdgraph.obj3 ) ) main = c( main, "estimate3" )
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

		row_plot = ifelse( is.null( bdgraph.obj2 ), 1, 2 )
		op       = graphics::par( mfrow = c( row_plot, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )

		igraph::plot.igraph( G_igraph,   layout = igraph::layout.circle, main = main[1], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		igraph::plot.igraph( est_igraph, layout = igraph::layout.circle, main = main[2], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		 
		if( !is.null( bdgraph.obj2 ) )
		{
			est2_igraph <- igraph::graph.adjacency( as.matrix( est2 ), mode = "undirected", diag = FALSE )
			igraph::plot.igraph( est2_igraph, layout = igraph::layout.circle, main = main[3], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
		}
		
		if( !is.null( bdgraph.obj3 ) )
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
		
	tp = sum( ( upper_G != 0 ) * ( upper_est_G != 0 ) ) 
	tn = sum( ( upper_G == 0 ) * ( upper_est_G == 0 ) )
	fp = sum( ( upper_G == 0 ) * ( upper_est_G != 0 ) ) 
	fn = sum( ( upper_G != 0 ) * ( upper_est_G == 0 ) )
		
	# harmonic mean of precision and recall, called F-measure or balanced F-score
	F1score = ( 2 * tp ) / ( 2 * tp + fp + fn )

	specificity  = tn / ( tn + fp )
	sensitivity  = tp / ( tp + fn )
	# Matthews Correlation Coefficients (MCC)
	mcc          = ( ( tp * tn ) - ( fp * fn ) ) / ( sqrt( ( tp + fp ) * ( tp + fn ) ) * sqrt( ( tn + fp ) * ( tn + fn ) ) )
	
	return( c( tp, tn, fp, fn, F1score, specificity, sensitivity, mcc ) )
}
   
