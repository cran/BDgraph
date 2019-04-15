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
    G   = BDgraph::get_graph( target )
    est = BDgraph::get_graph( est    )
    
    p = ncol( G )
    result                 = matrix( 1, 8, 2 )
    result[ 1        , 1 ] = sum( G[ upper.tri( G ) == 1 ] )
    result[ 2        , 1 ] = p * ( p - 1 ) / 2 - result[ 1, 1 ]  
    result[ c( 3, 4 ), 1 ] = 0
    
    result[ , 2 ] = compute_measures( G = G, est_G = est )
    
    if( !is.null( est2 ) )
    {
        est2 = BDgraph::get_graph( est2 )
        result = cbind( result, compute_measures( G = G, est_G = est2 ) )
    }
    
    if( !is.null( est3 ) )
    {
        est3 = BDgraph::get_graph( est3 )
        result = cbind( result, compute_measures( G = G, est_G = est3 ) )
    }
    
    if( !is.null( est4 ) )
    {
        est4 = BDgraph::get_graph( est4 )
        result = cbind( result, compute_measures( G = G, est_G = est4 ) )
    }
    
    result[ is.na( result ) ] = 0
    
    if( is.null( main ) ) 
    {
        main = c( "Target", "estimate1" )
        if( !is.null( est2 ) ) main = c( main, "estimate2" )
        if( !is.null( est3 ) ) main = c( main, "estimate3" )
        if( !is.null( est4 ) ) main = c( main, "estimate4" )
    }
    
    colnames( result ) <- main
    rownames( result ) <- c( "true positive", "true negative", "false positive", "false negative", 
                             "F1-score", "specificity", "sensitivity", "MCC" )
    
    if( vis == TRUE )
    {
        row_plot = ifelse( is.null( est2 ), 1, 2 )
        op       = graphics::par( mfrow = c( row_plot, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )
        
        BDgraph::plot.graph( G  , main = main[1] )
        BDgraph::plot.graph( est, main = main[2] )

        if( !is.null( est2 ) ) BDgraph::plot.graph( est2, main = main[3] )
        if( !is.null( est3 ) ) BDgraph::plot.graph( est3, main = main[4] )
        if( !is.null( est4 ) ) BDgraph::plot.graph( est4, main = main[5] )

        graphics::par( op )
    }
    
    return( round( result, 3 ) )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    To compare measures the performance of estimated graphs based on true graph
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
compute_measures = function( G, est_G ) 
{
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
