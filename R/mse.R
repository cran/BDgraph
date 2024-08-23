## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2023  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Computes (weighted) mean squared error.                                  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

mse = function( pred, actual, weight = FALSE )
{
    if( ( inherits( pred, "bdgraph" ) ) || ( inherits( pred, "ssgraph" ) ) )
        pred = plinks( pred, round = 10 )
    
    if( !is.vector( actual ) )  
        actual = BDgraph::get_graph( actual )
    
    if( is.matrix( actual ) )
        actual_vec = actual[ upper.tri( actual ) ]
    
    if( is.matrix( pred ) )
        pred_vec = pred[ upper.tri( pred ) ]
    
    if( is.numeric( weight ) )
    {
        actual_one  = actual_vec[ actual_vec == 1 ]
        actual_zero = actual_vec[ actual_vec == 0 ]
        
        pred_one  = pred_vec[ actual_vec == 1 ]
        pred_zero = pred_vec[ actual_vec == 0 ]
      
        mse_one  = mean( ( actual_one - pred_one ) ^ 2 )
        mse_zero = mean( ( actual_zero - pred_zero ) ^ 2 )
        
        mse_value = weight * mse_one +  ( 1 - weight ) * mse_zero 
         
    }else{
      
        if( weight == FALSE )
            mse_value = mean( ( actual_vec - pred_vec ) ^ 2 )
    }
    
    return( mse_value )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
