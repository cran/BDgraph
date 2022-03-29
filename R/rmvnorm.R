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
#     Data generator from multivarate normal distribution                      |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

rmvnorm = function( n = 10, mean = rep( 0, length = ncol( sigma ) ), 
                    sigma = diag( length( mean ) ) )
{
    if( !isSymmetric( sigma, tol = sqrt( .Machine $ double.eps ), check.attributes = FALSE ) ) 
        stop( "'sigma' must be a symmetric matrix" )
    
    sigma = as.matrix( sigma )
    p     = nrow( sigma )
    
    if( length( mean ) == 1 ) mean = rep( mean, p )
    if( length( mean ) != nrow( sigma ) ) stop( "'mean' and 'sigma' have non-conforming size" )
    
    # - - generate multivariate normal data - - - - - - - - - - - - - - - - - -|
    
    chol_sig = chol( sigma )
    z        = matrix( stats::rnorm( p * n ), p, n )
    data     = t( chol_sig ) %*% z + mean
    data     = t( data )
    
    if( n == 1 ) data = as.vector( data )
        
    return( data )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
