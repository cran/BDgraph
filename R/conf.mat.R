## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2021  Reza Mohammadi                                |
#                                                                              |
#     This file is part of 'BDgraph' package.                                  |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Create a Confusion Matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

conf.mat = function( pred, actual, cutoff = 0.5, proportion = FALSE, 
                     dnn = c( "Prediction", "Actual" ), ... )
{
    if( is.matrix( actual ) ) 
        if( ( sum( actual == 0 ) + sum( actual == 1 ) ) != ( nrow( actual ) ^ 2 ) ) stop( "Elements of matrix 'actual' must be 0 or 1" )
    
    if( inherits( actual, "sim" ) ) actual = actual $ G
    
    if( is.matrix( pred ) )
        if( any( pred < 0 ) || any( pred > 1 ) ) stop( "Elements of matrix 'pred' must be between ( 0, 1 )" )
    
    if( ( inherits( pred, "bdgraph" ) ) | ( inherits( pred, "ssgraph" ) ) )
    {
        pred = pred $ p_links
        if( is.null( pred ) ) pred = BDgraph::plinks( pred, round = 15 )
    }
    
    pred   = pred[   upper.tri( pred   ) ]    
    actual = actual[ upper.tri( actual ) ]   
    
    if( length( pred ) != length( actual ) ) stop( "'prod' and 'actual' lengths differ" )
    
    if( ( cutoff < 0 ) || ( cutoff > 1 ) ) stop( "'cutoff' must be between 0 and 1" )
        
    pred = ifelse( pred >= cutoff, 1, 0 )
        
    conf_mat = table( pred, actual, dnn = dnn, ... )
    
    if( proportion == TRUE )
        conf_mat = round( conf_mat / sum( conf_mat ), 3 )
    
    return( conf_mat )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |







