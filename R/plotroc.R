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
#     To plot ROC curve                                                        |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

plotroc = function( pred, actual, cut = 200, smooth = FALSE, calibrate = TRUE, 
                    linetype = NULL, color = NULL, size = 1, main = "ROC Curve", 
                    xlab = "False Postive Rate", ylab = "True Postive Rate",
                    legend = TRUE, legend.size = 17, legend.position = c( 0.7, 0.3 ),
                    labels = NULL, auc = TRUE, theme = ggplot2::theme_minimal() )
{
    if( !inherits( pred, "list" ) ) 
        pred = list( pred )
    
    length_pred = length( pred )
    
    if( is.null( color    ) ) color    = 1:length_pred
    if( is.null( linetype ) ) linetype = 1:length_pred
    
    if( !is.vector( actual ) )
    {
        adj_G  = BDgraph::get_graph( actual )
        
        actual = adj_G[ upper.tri( adj_G ) ]
    }
    
    if( legend && is.null( labels ) )
    {
        labels = numeric( length = length_pred )
        
        for( i in 1:length_pred ) 
            labels[ i ] = paste0( "pred ", i )
    }
    
    fp_roc = vector()
    tp_roc = vector()
    
    length_fp_i = numeric( length = length_pred )
    
    for( i in 1:length_pred )
    {
        output_tp_fp = compute_tp_fp( pred = pred[[i]], actual = actual, cut = cut, smooth = smooth, calibrate = calibrate )
        
        length_fp_i[ i ] = length( output_tp_fp $ fp )
        
        fp_roc = c( fp_roc, sort( c( output_tp_fp $ fp ) ) )
        tp_roc = c( tp_roc, sort( c( output_tp_fp $ tp ) ) )
    }

    if( legend && auc )
    {
        auc_s = BDgraph::auc( pred = pred, actual = actual, cut = cut, calibrate = calibrate )
            
        labels = paste( labels, "; AUC=", round( auc_s, 3 ) )
    }
      
    df_gg = data.frame( pred = rep( as.factor( 1:length_pred ), length_fp_i ), 
                        fp_roc = fp_roc, tp_roc = tp_roc )

    ggplot2::ggplot( df_gg, ggplot2::aes( x = fp_roc, y = tp_roc ) ) +
        ggplot2::geom_line( ggplot2::aes( group = pred, colour = pred, linetype = pred ), show.legend = legend, size = size ) +
        ggplot2::scale_color_manual( values = color, labels = labels ) +
        ggplot2::scale_linetype_manual( values = linetype, labels = labels ) +
        ggplot2::labs( x = xlab, y = ylab ) +
        ggplot2::ggtitle( main ) +
        theme +
        ggplot2::theme( legend.title = ggplot2::element_blank(), legend.position = legend.position, 
               text = ggplot2::element_text( size = legend.size ) ) +
        ggplot2::theme( legend.key.width = ggplot2::unit( 2, "line" ) ) 
}
       
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Function to compute tp (true positive) and fp (false positive) for ROC plot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

compute_tp_fp = function( pred, actual, cut, smooth, calibrate = TRUE )
{
    if( !is.vector( actual ) )
    {
        adj_G  = BDgraph::get_graph( actual )
        actual = adj_G[ upper.tri( adj_G ) ]
    }

    if( ( inherits( pred, "bdgraph" ) ) | ( inherits( pred, "ssgraph" ) ) )
        pred = BDgraph::plinks( pred, round = 15 )
    
    if( is.matrix( pred ) )
    {
        if( any( round( pred ) < 0 ) || any( round( pred ) > 1 ) ) 
            stop( "Elements of 'pred' must be between [ 0, 1 ]" )
        
        pred = pred[ upper.tri( pred ) ]
    }

    sum_one  = sum( actual )
    sum_zero = length( actual ) - sum_one
        
    if( !inherits( pred, "huge" ) )
    {
        tp = c( 1, rep( 0, cut ) )
        fp = tp
        
        cut_points = ( 0 : cut ) / cut

        if( calibrate == FALSE ) w_roc = 1
        
        for( i in 2 : cut )
        {
            if( calibrate == TRUE ) w_roc = abs( pred - cut_points[ i ] )

            tp[ i ] = sum( w_roc * ( actual == 1 ) * ( pred > cut_points[ i ] ) ) / sum_one
            fp[ i ] = sum( w_roc * ( actual == 0 ) * ( pred > cut_points[ i ] ) ) / sum_zero
        }
    }
    
    if( inherits( pred, "huge" ) )
    {
        path = pred $ path
        tp   = numeric( length( path ) )
        fp   = tp
        
        for( i in 1 : length( path ) )
        {
            est_G  = as.matrix( path[[ i ]] )
            pred_i = est_G[ upper.tri( est_G ) ]
            
            tp[ i ] = sum( ( actual == 1 ) * ( pred_i == 1 ) ) / sum_one
            fp[ i ] = sum( ( actual == 0 ) * ( pred_i == 1 ) ) / sum_zero
        }
        
        tp = c( tp, 1 )
        fp = c( fp, 1 )
    }
    
    if ( smooth == TRUE )
    {
        fit = stats::smooth.spline( x = fp, y = tp )
        fp  = c( 0, fit $ x )
        tp  = c( 0, fit $ y )
    }	
    
    return( list( tp = tp, fp = fp ) )
}
  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
