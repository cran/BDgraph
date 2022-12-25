## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2022  Reza Mohammadi                                |
#                                                                              |
#     This file is part of BDgraph package.                                    |
#                                                                              |
#     BDgraph is free software: you can redistribute it and/or modify it under |
#     the terms of the GNU General Public License as published by the Free     |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.|
#                                                                              |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    posterior predict function for "bdgraph" object                           |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

posterior.predict = function( object, iter = 1, ... )
{
    if( is.null( object $ all_graphs ) ) 
        stop( "'bdgraph.obj' must be an object of function 'bdgraph()' or 'ssgraph()' with option 'save = TRUE'" )
    
    method = object $ method
    data   = object $ data
    
    n_data = nrow( data )
    p      = ncol( data )

	if( isSymmetric( data ) )
	{
		S = data
	}else{
 		S = t( data ) %*% data
	}

	sample_graphs = object $ sample_graphs
    all_graphs    = object $ all_graphs
	graph_weights = object $ graph_weights
        
    sample_G = sample( x = sample_graphs, size = iter, replace = TRUE, prob = graph_weights ) 
    
    G_i = matrix( 0, nrow = p, ncol = p )
    upper_G_i = G_i[ upper.tri( G_i ) ]
    
    Z = matrix( 0, nrow = iter, ncol = p )
        
    for( i in 1:iter )
    {
        upper_G_i = upper_G_i * 0
        
        upper_G_i[ which( unlist( strsplit( as.character( sample_G[i] ), "" ) ) == 1 ) ] = 1
        
        G_i[ upper.tri( G_i ) ] = upper_G_i

        K_i = BDgraph::rgwish( n = 1, adj = G_i, b = 3 + n_data, D = diag( p ) + S )
        sigma_i = solve( K_i )
        
        Z[ i, ] = BDgraph::rmvnorm( n = 1, mean = 0, sigma = sigma_i )
    }
    
    if( method == "ggm" )
        sample = Z

	if( method == "tgm" )
	{
	    mean = 0
	    nu   = 1
	    
	    tau_gamma = stats::rgamma( n = iter, shape = nu / 2, rate = nu / 2 )
        sample    = mean + Z / sqrt( tau_gamma )
	}
    
    if( method == "gcgm" ) 
    {
        K = object $ K_hat
        
        if( is.null( K ) )
        {
            G = BDgraph::select( bdgraph.obj = object )
            
            sample_K = BDgraph::rgwish( n = 500, adj = G, b = 3 + n_data, D = diag( p ) + S )
            
            K = 0 * G
            for( i in 1:dim( sample_K )[3] )
                K = K + sample_K[[i]]
            
            K = K / dim( sample_K )[3]
        }
        
        sample = 0 * Z
        
        for( j in 1:p )
        {
            sdj = sqrt( 1 / K[ j, j ] )     # 2a: # variance of component j (given the rest!)
            muj = - sum( Z[ , -j, drop = FALSE ] %*% K[ -j, j, drop = FALSE ] / K[ j, j ] )	 
            
            table_j = table( data[ , j ] )
            cat_y_j = as.numeric( names( table_j ) ) 
            len_cat_y_j = length( cat_y_j )
            
            if( len_cat_y_j > 1 )
            {
                cum_prop_yj = cumsum( table_j[ -len_cat_y_j ] ) / n_data
                
                #cut_j = vector( length = len_cat_y_j - 1 )
                # for( k in 1:length( cut_j ) ) cut_j[ k ] = stats::qnorm( cum_prop_yj[ k ] )
                cut_j = stats::qnorm( cum_prop_yj, mean = 0, sd = 1 )
                            
            	breaks = c( min( Z[ , j ] ) - 1, cut_j, max( Z[ , j ] ) + 1 )  
            	
            	ind_sj = as.integer( cut( Z[ , j ], breaks = breaks, right = FALSE ) )
            	
            	sample[ , j ]  = cat_y_j[ ind_sj ]
            }else{
                sample[ , j ]  = cat_y_j
            }
        }
    }

    if( method == "dw" )
    {
        q    = object $ q.est
        beta = object $ beta.est
        mean = rep( 0, p )
        
        #Z = tmvtnorm::rtmvnorm( n = iter, mean = mean, sigma = sigma, lower = rep( -5, length = p ), upper = rep( 5, length = p ) )
        
        pnorm_Z = stats::pnorm( Z )
        
        if( is.matrix( q ) && is.matrix( beta ) )
        {
            for( j in 1 : p ) 
                sample[ ,j ] = BDgraph::qdweibull( pnorm_Z[ , j ], q = q[ , j ], beta = beta[ , j ], zero = TRUE )
        }
        
        if( is.vector( q ) && is.vector( beta ) )
        {
            for( j in 1 : p ) 
                sample[ , j ] = BDgraph::qdweibull( pnorm_Z[ , j ], q = q[ j ], beta = beta[ j ], zero = TRUE )		    
        }
    }
        
    return( sample )
}
  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |














