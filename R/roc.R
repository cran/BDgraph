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

roc = function( pred, actual, auc = TRUE, smooth = FALSE, plot = FALSE, quiet = TRUE, ... )
{
    if( is.matrix( actual ) ) 
        if( ( sum( actual == 0 ) + sum( actual == 1 ) ) != ( nrow( actual ) ^ 2 ) ) stop( "Elements of matrix 'actual' must be 0 or 1" )
    
    if( inherits( actual, "sim" ) ) actual = actual $ G
    
    response = actual[ upper.tri( actual ) ]   
    
    if( is.matrix( pred ) )
        if( any( pred < 0 ) || any( pred > 1 ) ) stop( "Elements of matrix 'pred' must be between ( 0, 1 )" )
    
    if( ( inherits( pred, "bdgraph" ) ) | ( inherits( pred, "ssgraph" ) ) )
        pred = BDgraph::plinks( pred, round = 15 )
    
    predictor = pred[ upper.tri( pred ) ]    
    
    pROC::roc( response = response, predictor = predictor, levels = c( 0, 1 ), 
               quiet = quiet, smooth = smooth, plot = plot, ... )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |




