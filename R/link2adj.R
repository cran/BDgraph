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
#     Creating an adjacency matrix based on links                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

link2adj = function( link, p = NULL )
{
    if( !is.matrix( link ) & !is.data.frame( link ) ) stop( "Input 'link' must be a matrix or dataframe" )
    if( is.data.frame( link ) ) link <- data.matrix( link )
    
    if( ncol( link ) != 2 ) stop( "'link' must have only 2 columns" )
    if( nrow( link ) < 1  ) stop( "'link' must have at least one row" )
    
    if( !is.null( p ) ) 
        if( max( link ) > p ) stop( "'p' is not matched with input 'link'" )
    
    if( is.null( p ) ) p = max( link )
    
    adj = matrix( 0, p, p )
    
    for( i in 1:nrow( link ) )
        adj[ link[ i, 1 ], link[ i, 2 ] ] = 1 
    
    return( adj ) 
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
