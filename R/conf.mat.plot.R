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
#     Plot a Confusion Matrix
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

conf.mat.plot = function( pred, actual, cutoff = 0.5, conf.level = 0, margin = 1, 
                          color = c( "#ff83a8", "#83ff9b" ), ... )
{
    conf_mat = BDgraph::conf.mat( pred = pred, actual = actual, cutoff = cutoff )
    
    graphics::fourfoldplot( conf_mat, conf.level = conf.level, 
                            margin = margin, color = color, ... )
}   
  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |




