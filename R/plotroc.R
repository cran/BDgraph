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

plotroc = function( actual, pred, cut = 20, smooth = FALSE, label = TRUE, 
                    AUC = TRUE, main = "ROC Curve" )
{
    if( !inherits( pred, "list" ) ) pred = list( pred )
    length_pred = length( pred )
  
    if( length_pred == 1 ) label = FALSE
    
    G = BDgraph::get_graph( actual )  
    G[ lower.tri( G, diag = TRUE ) ] = 0
    
    if( length( label ) > 1 ){
        labels_plot = label
        label = TRUE
    }else{
        labels_plot = rep( NA, length_pred )
        if( label != TRUE ) label = FALSE
    }
    
    fp_roc = vector()
    tp_roc = vector()
    
    length_fp_i = numeric( length = length_pred )
    
    for( i in 1:length_pred )
    {
        output_tp_fp = compute_tp_fp( G = G, est = pred[[i]], cut = cut, smooth = smooth )
        
        length_fp_i[ i ] = length( output_tp_fp $ fp )
        
        fp_roc = c( fp_roc, sort( c( output_tp_fp $ fp ) ) )
        tp_roc = c( tp_roc, sort( c( output_tp_fp $ tp ) ) )
        
        if( label && is.na( labels_plot[ i ] ) ) labels_plot[ i ] = paste0( "pred ", i )
    }

    if( label && AUC )
    {
        for( i in 1:length_pred ) 
        {
            roc_i = BDgraph::roc( actual = G, pred = pred[[i]] )
            
            labels_plot[ i ] = paste( labels_plot[ i ], "; AUC=", round( pROC::auc( roc_i ), 3 ) )
        }
    }
      
    df_gg = data.frame( pred = rep( as.factor( 1:length_pred ), length_fp_i ), 
                        fp_roc = fp_roc, tp_roc = tp_roc )

    ggplot2::ggplot( df_gg, ggplot2::aes( x = fp_roc, y = tp_roc ) ) +
        ggplot2::geom_line( ggplot2::aes( group = pred, colour = pred, linetype = pred ), show.legend = label ) +
        ggplot2::scale_color_manual( values = 1:length_pred, labels = labels_plot ) +
        ggplot2::scale_linetype_manual( values = 1:length_pred, labels = labels_plot ) +
        ggplot2::labs( x = "False Postive Rate", y = "True Postive Rate" ) +
        ggplot2::ggtitle( main ) +
        ggplot2::theme_minimal() +
        ggplot2::theme( legend.title = ggplot2::element_blank(), legend.position = c( .7, .3 ), 
               text = ggplot2::element_text( size = 17 ) ) +
        ggplot2::theme( legend.key.width = ggplot2::unit( 2, "line" ) ) 
}
       
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Function to compute tp (true positive) and fp (false positive) for ROC plot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
compute_tp_fp = function( G, est, cut, smooth )
{
    p           = nrow( G )
    upper_G     = G[ upper.tri( G ) ]
    sum_edges   = sum( upper_G )
    sum_no_dges = p * ( p - 1 ) / 2 - sum_edges
    
    if( ( inherits( est, "bdgraph" ) ) | ( inherits( est, "ssgraph" ) ) )
    {
        p_links = est $ p_links
        if( is.null( p_links ) ) p_links = BDgraph::plinks( est, round = 15 )
    }
    
    if( is.matrix( est ) )
    {
        if( any( est < 0 ) || any( est > 1 ) ) stop( "Elements of 'est' must be between ( 0, 1 )" )
        p_links = est
    }
    
    if( !inherits( est, "huge" ) )
    {
        tp = c( 1, rep( 0, cut ) )
        fp = tp
        
        cut_points = ( 0 : cut ) / cut
        
        for( i in 2 : cut )
        {
            # checking for cut pints
            est_G = matrix( 0, p, p )
            est_G[ p_links > cut_points[ i ] ] = 1
            upper_est_G = est_G[ upper.tri( est_G ) ]
            
            tp[ i ] = sum( ( upper_G == 1 ) * ( upper_est_G == 1 ) ) / sum_edges
            fp[ i ] = sum( ( upper_G == 0 ) * ( upper_est_G == 1 ) ) / sum_no_dges
        }
    }
    
    if( inherits( est, "huge" ) )
    {
        path = est $ path
        tp   = numeric( length( path ) )
        fp   = tp
        
        for( i in 1 : length( path ) )
        {
            est_G       = as.matrix( path[[ i ]] )
            upper_est_G = est_G[ upper.tri( est_G ) ]
            
            tp[ i ] = sum( ( upper_G == 1 ) * ( upper_est_G == 1 ) ) / sum_edges
            fp[ i ] = sum( ( upper_G == 0 ) * ( upper_est_G == 1 ) ) / sum_no_dges
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
