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

roc = function( pred, actual, auc = TRUE, smooth = FALSE, plot = FALSE, 
                quiet = TRUE, ... )
{
    if( ( inherits( pred, "bdgraph" ) ) | ( inherits( pred, "ssgraph" ) ) )
        pred = BDgraph::plinks( pred, round = 15 )

    if( inherits( actual, "sim" ) ) 
        actual = actual $ G
    
    if( is.matrix( pred   ) ) pred   = pred[   upper.tri( pred   ) ]
    if( is.matrix( actual ) ) actual = actual[ upper.tri( actual ) ]
    
    if( any( round( pred ) < 0 ) || any( round( pred ) > 1 ) ) 
        stop( "Elements of 'pred' must be between ( 0, 1 )" )
    
    if( ( sum( actual == 0 ) + sum( actual == 1 ) ) != length( actual ) ) 
        stop( "Elements of 'actual' must be 0 or 1" )
    
    pROC::roc( response = actual, predictor = pred, levels = c( 0, 1 ), 
               quiet = quiet, smooth = smooth, plot = plot, ... )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |




