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
#     The Discrete Weibull Distribution ( Type 1 )                             |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

ddweibull = function( x, q = exp( -1 ), beta = 1, zero = TRUE )
{
    if( any( x != floor( x ) )      ) stop( "'x' must be an integer" )
    if( max( q ) > 1 | min( q ) < 0 ) stop( "'q' must be between 0 and 1" )
    if( min( beta ) <= 0            ) stop( "'beta' must be a positive value" )
    
    if( zero == FALSE ) x = x - 1
    
    return( q ^ x ^ beta - q ^ ( x + 1 ) ^ beta )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
pdweibull = function( x, q = exp( -1 ), beta = 1, zero = TRUE )
{
    if( max( q ) > 1 | min( q ) < 0 ) stop( "'q' must be between 0 and 1" )
    if( min( beta ) <= 0 ) stop( "'beta' must be a positive value" )
    
    if( zero == FALSE ) x = x - 1
    
    p_dw = 1 - q ^ ( x + 1 ) ^ beta
    p_dw[ x < 0 ] = 0
    
    return( p_dw )
 } 
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
qdweibull = function( p, q = exp( -1 ), beta = 1, zero = TRUE )
{
    if( max( q ) > 1 | min( q ) < 0 ) stop( "'q' must be between 0 and 1" )
    if( min( beta ) <= 0 ) stop( "'beta' must be a positive value" )
    
    if( zero )
        return( ceiling( ( log( 1 - p ) / log( q ) ) ^ ( 1 / beta ) - 1 ) )
    else
        return( ceiling( ( log( 1 - p ) / log( q ) ) ^ ( 1 / beta ) ) )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
rdweibull = function( n, q = exp( -1 ), beta = 1, zero = TRUE )
{
    if( max( q ) > 1 | min( q ) < 0 ) stop( "'q' must be between 0 and 1" )
    if( min( beta ) <= 0 ) stop( "'beta' must be a positive value" )
    
    r_unif = stats::runif( n )
    
    if( zero )
        return( ceiling( ( log( 1 - r_unif ) / log( q ) ) ^ ( 1 / beta ) ) - 1 )
    else
        return( ceiling( ( log( 1 - r_unif ) / log( q ) ) ^ ( 1 / beta ) ) )
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |






