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
#     Reports the measures to assess the performance of estimated graphs       |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

compare = function( pred, actual, main = NULL, vis = FALSE ) 
{
    if( !inherits( pred, "list" ) ) pred = list( pred )
    length_pred = length( pred )
    
    G = BDgraph::get_graph( actual )
    
    p = ncol( G )
    
    result                 = matrix( 1, 8, length_pred + 1 )
    result[ 1        , 1 ] = sum( G[ upper.tri( G ) == 1 ] )
    result[ 2        , 1 ] = p * ( p - 1 ) / 2 - result[ 1, 1 ]  
    result[ c( 3, 4 ), 1 ] = 0
    
    for( i in 1:length_pred )
    {
        est_g = BDgraph::get_graph( pred[[i]] )
        
        result[ , i + 1 ] = compute_measures( est_G = est_g, G = G )
    }
    
    result[ is.na( result ) ] = 0
    
    rownames( result ) <- c( "True Positive", "True Negative", "False Positive", "False Negative", 
                             "F1-score", "Specificity", "Sensitivity", "MCC" )
    
    if( is.null( main ) ) 
    {
        main = c( "Actual" )
        
        for( i in 1:length_pred )
            main = c( main, paste0( "pred ", i ) )
    }
    
    colnames( result ) = main

    if( vis == TRUE )
    {
        if( length_pred == 1 ) mfrow = c( 1, 2 )  
        if( ( length_pred > 1  ) & ( length_pred < 4  ) ) mfrow = c( 2, 2 )  
        if( ( length_pred > 3  ) & ( length_pred < 6  ) ) mfrow = c( 3, 2 )  
        if( ( length_pred > 5  ) & ( length_pred < 9  ) ) mfrow = c( 3, 3 )  
        if( ( length_pred > 8  ) & ( length_pred < 12 ) ) mfrow = c( 4, 3 )  
        if( ( length_pred > 11 ) & ( length_pred < 16 ) ) mfrow = c( 4, 4 )  
        
        op = graphics::par( mfrow = mfrow, pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) )
        
        BDgraph::plot.graph( G  , main = main[ 1 ] )
        
        for( i in 1:length_pred )
            BDgraph::plot.graph( pred[[i]], main = main[ i + 1 ] )

        graphics::par( op )
    }
    
    return( round( result, 3 ) )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    To compare measures the performance of estimated graphs based on true graph
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
compute_measures = function( est_G, G ) 
{
    if( sum( dim( G ) == dim( est_G ) ) != 2 ) 
        stop( "'pred' and 'actual' have non-conforming size" )
    
    upper_G     = G[     upper.tri( G     ) ]
    upper_est_G = est_G[ upper.tri( est_G ) ]
    
    tp = sum( ( upper_G == 1 ) * ( upper_est_G == 1 ) )  # True  Positive
    tn = sum( ( upper_G == 0 ) * ( upper_est_G == 0 ) )  # True  Negative
    fp = sum( ( upper_G == 0 ) * ( upper_est_G == 1 ) )  # False Positive
    fn = sum( ( upper_G == 1 ) * ( upper_est_G == 0 ) )  # False Negative
    
    # harmonic mean of precision and recall, called F-measure or balanced F-score
    F1score = ( 2 * tp ) / ( 2 * tp + fp + fn )
    
    specificity = tn / ( tn + fp )
    sensitivity = tp / ( tp + fn )
    
    # Matthews Correlation Coefficients (MCC)
    mcc = ( ( tp * tn ) - ( fp * fn ) ) / ( sqrt( ( tp + fp ) * ( tp + fn ) ) * sqrt( ( tn + fp ) * ( tn + fn ) ) )
    
    return( c( tp, tn, fp, fn, F1score, specificity, sensitivity, mcc ) )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
