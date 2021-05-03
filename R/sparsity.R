## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2020  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Compute the sparsity of an adjacency matrix                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

sparsity = function( adj )
{
    G = BDgraph::get_graph( adj )
    p = ncol( G )
    
    sum_E = sum( adj ) / 2
    D     = p * ( p - 1 ) / 2
    
    sparsity_g = sum_E / D
    
    return( sparsity_g )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
