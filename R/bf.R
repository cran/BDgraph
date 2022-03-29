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
#     Computing the Bayes factor between two graph structures                  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bf = function( num, den, bdgraph.obj, log = TRUE )
{
    G_1 = BDgraph::get_graph( num )
    G_2 = BDgraph::get_graph( den )
    
    if( ( inherits( bdgraph.obj, "bdgraph" ) ) | ( inherits( bdgraph.obj, "ssgraph" ) ) )
    {
        p_links = bdgraph.obj $ p_links
        if( is.null( p_links ) ) p_links = BDgraph::plinks( bdgraph.obj, round = 15 )
        min_plinks = 1 - p_links
        
        p_links[    p_links    == 0 ] = .Machine $ double.xmin
        min_plinks[ min_plinks == 0 ] = .Machine $ double.xmin
        
        g1_min_g2 = G_1 - G_2
        
        upper_plinks     = p_links[    upper.tri( p_links    ) ]
        upper_min_plinks = min_plinks[ upper.tri( min_plinks ) ]
        
        upper_g1_min_g2 = g1_min_g2[ upper.tri( g1_min_g2 ) ]

        log_vec_bf = 0 * upper_plinks
        log_vec_bf[ upper_g1_min_g2 == 1  ] = upper_plinks[      upper_g1_min_g2 == 1  ] / upper_min_plinks[ upper_g1_min_g2 == 1  ]
        log_vec_bf[ upper_g1_min_g2 == -1 ] =  upper_min_plinks[ upper_g1_min_g2 == -1 ] / upper_plinks[     upper_g1_min_g2 == -1 ]
        log_vec_bf = log_vec_bf[ log_vec_bf != 0 ]
        
        if( log == TRUE )
        {
            bf = sum( log( log_vec_bf ) )
        }else{
            bf = prod( log_vec_bf )
        }
    }
    
    return( bf )
}
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |






