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
#     Extract links from an adjacency matrix                                   |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

adj2link = function( adj )
{
    G = BDgraph::get_graph( adj )
    if( sum( G ) == 0 ) print( " 'adj' has no link. " )
    p = ncol( G )
    
    links = vector()
    for( i in 1:p )
        for( j in 1:p )
            if( G[ i, j ] == 1 ) links = rbind( links, c( i, j ) )
    
    return( links ) 
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
