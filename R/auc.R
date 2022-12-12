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

auc = function( pred, actual, cut = 200, calibrate = TRUE )
{
    if( !inherits( pred, "list" ) ) 
        pred = list( pred )
    
    length_pred = length( pred )
    
    if( !is.vector( actual ) )
    {
        adj_G  = BDgraph::get_graph( actual )
        
        actual = adj_G[ upper.tri( adj_G ) ]
    }
  
    if( ( sum( actual == 0 ) + sum( actual == 1 ) ) != length( actual ) ) 
        stop( "Elements of 'actual' must be 0 or 1" )
    
    auc_value = vector()
    
    for( i in 1:length_pred )
    {
        output_tp_fp = compute_tp_fp( pred = pred[[i]], actual = actual, cut = cut, smooth = FALSE, calibrate = calibrate )
            
        fp = output_tp_fp $ fp
        tp = output_tp_fp $ tp
            
        diffs_x    =  fp[ -length( fp ) ] - fp[ -1 ]
        means_vert = ( tp[ -1 ] + tp[ -length( tp ) ] ) / 2
        
        auc_value = c( auc_value, sum( means_vert * diffs_x ) )
    }
    
    return( auc_value )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |




