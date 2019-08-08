## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2019  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of BDgraph package.                                                        |
#                                                                                                  |
#     BDgraph is free software: you can redistribute it and/or modify it under                     |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Computing estimated percision matrix                                                                  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
covariance = function( bdgraph.obj, round = 2 )
{
    if( ( class( bdgraph.obj ) != "bdgraph" ) && ( class( bdgraph.obj ) != "ssgraph" ) )
        stop( "Input 'bdgraph.obj' must be an object from functions 'bdgraph()', 'bdgraph.mpl()', or 'ssgraph()'." )
    
    K_hat = bdgraph.obj $ K_hat
    
    if( is.null( K_hat ) )			  
        stop( " Input 'bdgraph.obj' must be an object from functions 'bdgraph()' or 'ssgraph()'." )
    
    cov = solve( K_hat )
    
    return( round( cov, round ) )
    
}  

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |













