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

compare = function( target, est, est2 = NULL, est3 = NULL, est4 = NULL, main = NULL, vis = FALSE ) 
{
    if( is.matrix( target ) ) 
    {
        if( ( sum( target == 0 ) + sum( target == 1 ) ) != ( nrow( target ) ^ 2 ) ) stop( "Element of 'target' must be 0 or 1" )
        G = target
    }
    
    if(   class( target ) == "sim"     ) G <- unclass( target $ G ) 
    if(   class( target ) == "graph"   ) G <- unclass( target ) 
    if( ( class( target ) == "bdgraph" ) | ( class( target ) == "ssgraph" ) ) G <- BDgraph::select( target ) 
    
    result        = matrix( 1, 8, 2 )
    result[ , 2 ] = compute_measures( G = G, est_G = est )
    
    if( !is.null( est2 ) )
        result = cbind( result, compute_measures( G = G, est_G = est2 ) )
    
    if( !is.null( est3 ) )
        result = cbind( result, compute_measures( G = G, est_G = est3 ) )
    
    if( !is.null( est4 ) )
        result = cbind( result, compute_measures( G = G, est_G = est4 ) )
    
    p = ncol( G )
    result[ c( 3, 4 ), 1    ] = 0
    result[ 1        , 1    ] = sum( G )
    result[ 2        , 1    ] = p * ( p - 1 ) / 2 - result[ 1, 1 ]  
    result[ is.na( result ) ] = 0
    
    if( is.null( main ) ) 
    {
        if( ( class( target ) == "bdgraph" ) | ( class( target ) == "ssgraph" ) )
        {
            main = c( "estimate1", "estimate2" )
            if( !is.null( est2 ) ) main = c( main, "estimate3" )
            if( !is.null( est3 ) ) main = c( main, "estimate4" )
            if( !is.null( est4 ) ) main = c( main, "estimate5" )
        }else{
            main = c( "True", "estimate" )
            if( !is.null( est2 ) ) main = c( main, "estimate2" )
            if( !is.null( est3 ) ) main = c( main, "estimate3" )
            if( !is.null( est4 ) ) main = c( main, "estimate4" )
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
        
        if( !is.null( est4 ) )
        { 
            est4_igraph <- igraph::graph.adjacency( as.matrix( est4 ), mode = "undirected", diag = FALSE ) 
            igraph::plot.igraph( est4_igraph, layout = igraph::layout.circle, main = main[5], vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )			
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
    if( is.matrix( est_G ) ) 
    {
        if( ( sum( est_G == 0 ) + sum( est_G == 1 ) ) != ( nrow( est_G ) ^ 2 ) ) stop( "Element of 'est' must be 0 or 1" )
    }else{
        if( ( class( est_G ) == "bdgraph" ) | ( class( est_G ) == "ssgraph" ) ) est_G <- BDgraph::select( est_G ) 
        if( class( est_G ) == "select"  )  est_G <- est_G $ refit
        est_G = as.matrix( est_G )
    }
    
    if( sum( dim( G ) == dim( est_G ) ) != 2 ) stop( " 'target' and 'est' have non-conforming size" )
    
    upper_G     = G[     upper.tri( G     ) ]
    upper_est_G = est_G[ upper.tri( est_G ) ]
    
    tp = sum( ( upper_G == 1 ) * ( upper_est_G == 1 ) )  # True  Positive
    tn = sum( ( upper_G == 0 ) * ( upper_est_G == 0 ) )  # True  Negative
    fp = sum( ( upper_G == 0 ) * ( upper_est_G == 1 ) )  # False Positive
    fn = sum( ( upper_G == 1 ) * ( upper_est_G == 0 ) )  # False Negative
    
    # harmonic mean of precision and recall, called F-measure or balanced F-score
    F1score = ( 2 * tp ) / ( 2 * tp + fp + fn )
    
    specificity  = tn / ( tn + fp )
    sensitivity  = tp / ( tp + fn )
    # Matthews Correlation Coefficients (MCC)
    mcc          = ( ( tp * tn ) - ( fp * fn ) ) / ( sqrt( ( tp + fp ) * ( tp + fn ) ) * sqrt( ( tn + fp ) * ( tn + fn ) ) )
    
    return( c( tp, tn, fp, fn, F1score, specificity, sensitivity, mcc ) )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
