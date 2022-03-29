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

plotroc = function( target, est, est2 = NULL, est3 = NULL, est4 = NULL, 
                    cut = 20, smooth = FALSE, label = TRUE, main = "ROC Curve" )
{
    if( is.matrix( target ) ) 
    {
        if( ( sum( target == 0 ) + sum( target == 1 ) ) != ( nrow( target ) ^ 2 ) ) stop( "Elements of 'target' must be 0 or 1" )
        G = target
    }
    
    if( inherits( target, "sim" )   ) G <- unclass( target $ G ) 
    if( inherits( target, "graph" ) ) G <- unclass( target ) 
    
    G[ lower.tri( G, diag = TRUE ) ] = 0
    
    output_tp_fp = compute_tp_fp( G = G, est = est, cut = cut, smooth = smooth )
    fp           = output_tp_fp $ fp
    tp           = output_tp_fp $ tp
    
    # graphics::par( mar = c( 3.8, 4.2, 1.8, 1 ) )
    graphics::plot( NA, type = "l", col = "black", cex.lab = 1.3, cex.main = 2, cex.axis = 1.2,
                    main = main, xlab = "False Postive Rate", ylab = "True Postive Rate", 
                    ylim = c( 0, 1 ), xlim = c( 0, 1 ) )
    graphics::points( x = fp, y = tp, type = "l", col = 1, lty = 1, lw = 2 )
    
    if( ( length( label ) == 1 ) && ( label == TRUE ) ) label = c( "est" )
    
    if( !is.null( est2 ) )
    {
        output_tp_fp = compute_tp_fp( G = G, est = est2, cut = cut, smooth = smooth )
        fp_2         = output_tp_fp $ fp
        tp_2         = output_tp_fp $ tp
        
        graphics::points( x = fp_2, y = tp_2, type = "l", col = 2, lty = 2, lw = 2 )
        if( ( length( label ) == 1 ) && ( label != FALSE ) ) label = c( label, "est2" )
    }
    
    if( !is.null( est3 ) )
    {   
        output_tp_fp = compute_tp_fp( G = G, est = est3, cut = cut, smooth = smooth )
        fp_3         = output_tp_fp $ fp
        tp_3         = output_tp_fp $ tp
        
        graphics::points( x = fp_3, y = tp_3, type = "l", col = 3, lty = 3, lw = 2 )
        if( ( length( label ) == 2 ) ) label = c( label, "est3" )
    }
    
    if( !is.null( est4 ) )
    {   
        output_tp_fp = compute_tp_fp( G = G, est = est4, cut = cut, smooth = smooth )
        fp_4         = output_tp_fp $ fp
        tp_4         = output_tp_fp $ tp
        
        graphics::points( x = fp_4, y = tp_4, type = "l", col = 4, lty = 4, lw = 2 )
        if( ( length( label ) == 3 ) ) label = c( label, "est3" )
    }
    
    if( any( label != FALSE ) ) 
        graphics::legend( "bottomright", label, lty = 1:length( label ), col = 1:length( label ), lwd = c( 2, 2 ), cex = 1.5 )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# Function to compute tp (true possitive) and fp (false possitive) for ROC plot
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
