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
#     Computing estimated precision matrix                                     |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
precision = function( bdgraph.obj, round = 2 )
{
    if( ( !inherits( bdgraph.obj, "bdgraph" ) ) && ( !inherits( bdgraph.obj, "ssgraph" ) ) )
        stop( "'bdgraph.obj' must be an object from functions 'bdgraph()', 'bdgraph.mpl()', or 'ssgraph()'" )
    
    K_hat = bdgraph.obj $ K_hat
    
    if( is.null( K_hat ) )			  
        stop( "Input 'object' must be from functions 'bdgraph()' or 'ssgraph()'" )

    return( round( K_hat, round ) )

}  
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |














