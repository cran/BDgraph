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
#     BDMCMC algorithm for graphical models based on Discrete Weibull          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

get_bounds_dw = function( data, q, beta, pii, n, p, zero = TRUE )
{
    lower_bounds = matrix( 0, nrow = n, ncol = p )
    upper_bounds = matrix( 0, nrow = n, ncol = p )
    
    if( is.vector( beta ) )
    {
        for ( j in 1:p ) 
        {   
            for( r in sort( unique( data[ , j ] ) ) )
            {
                ir = ( 1:n )[ data[ , j ] == r & !is.na( data[ , j ] ) ]
                
                pdw_lb = BDgraph::pdweibull( r - 1, q = q[ j ], beta = beta[ j ], zero = zero )
                pdw_ub = BDgraph::pdweibull( r    , q = q[ j ], beta = beta[ j ], zero = zero )
                
                lower_bounds[ ir, j ] = stats::qnorm( ( 1 - pii[ j ] ) * ( r != 0 ) + pii[ j ] * pdw_lb )
                upper_bounds[ ir, j ] = stats::qnorm( ( 1 - pii[ j ] ) + pii[ j ] * pdw_ub )
            }
        }
    }
    
    if( is.matrix( beta ) )
    {
        for ( j in 1:p ) 
        {   
            for( r in sort( unique( data[ , j ] ) ) )
            {
                ir = ( 1:n )[ data[ , j ] == r & !is.na( data[ , j ] ) ]

                pdw_lb = BDgraph::pdweibull( r - 1, q = q[ ir, j ], beta = beta[ ir, j ], zero = zero )
                pdw_ub = BDgraph::pdweibull( r    , q = q[ ir, j ], beta = beta[ ir, j ], zero = zero )
                
                lower_bounds[ ir, j ] = stats::qnorm( ( 1 - pii[ j ] ) * ( r != 0 ) + pii[ j ] * pdw_lb )
                upper_bounds[ ir, j ] = stats::qnorm( ( 1 - pii[ j ] ) + pii[ j ] * pdw_ub )
            }
        }
    }

    lower_bounds [ lower_bounds == -Inf ] = - .Machine $ double.xmax
    upper_bounds [ upper_bounds ==  Inf ] =   .Machine $ double.xmax

    lower_bounds [ lower_bounds ==  Inf ] =   .Machine $ double.xmax
    upper_bounds [ upper_bounds == -Inf ] = - .Machine $ double.xmax
    
    return( list( lower_bounds = lower_bounds, upper_bounds = upper_bounds ) )
}
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bdgraph.dw = function( data, x = NULL, formula = y ~ ., 
                       n = NULL, algorithm = "bdmcmc", iter = 5000, 
                       burnin = iter / 2, g.prior = 0.2, df.prior = 3,
                       ZI = FALSE, iter_bdw = 5000,
                       g.start = "empty",jump = NULL, save = FALSE, 
                       q = NULL, beta = NULL, pii = NULL,
                       cores = NULL, threshold = 1e-8, verbose = TRUE )
{
    if( is.matrix( data ) | is.data.frame( data ) ) 
        if( any( data < 0 ) ) 
            stop( "'data' should not have negative values" )
        
    if( df.prior < 3  ) stop( "'prior.df' must be >= 3" )
    if( iter < burnin ) stop( "'iter' must be higher than 'burnin'" )
    
    burnin = floor( burnin )
    
    if( is.numeric( verbose ) )
    {
        if( ( verbose < 1 ) | ( verbose > 100 ) ) 
            stop( "'verbose' (for numeric case) must be between ( 1, 100 )" )
        
        trace_mcmc = floor( verbose )
        verbose = TRUE
    }else{
        trace_mcmc = ifelse( verbose == TRUE, 10, iter + 1000 )
    }

    list_S_n_p = BDgraph::get_S_n_p( data = data, method = "dw", n = n, not.cont = NULL )
    
    S      = list_S_n_p $ S
    n      = list_S_n_p $ n
    p      = list_S_n_p $ p
    method = list_S_n_p $ method
    colnames_data = list_S_n_p $ colnames_data

    if( ( is.null( cores ) ) & ( p < 16 ) ) 
        cours = 1
        
    cores = BDgraph::get_cores( cores = cores, verbose = verbose )

    not.cont = list_S_n_p $ not.cont
    Z        = list_S_n_p $ Z
    data     = list_S_n_p $ data
    gcgm_NA  = list_S_n_p $ gcgm_NA

    sample_marginals = NULL
    
    if( is.null( q    ) ) if( inherits( data, "sim" ) ) q    = data $ q
    if( is.null( beta ) ) if( inherits( data, "sim" ) ) beta = data $ beta
    
    if( is.null( pii ) ) pii = rep( 1, p )

    if( is.null( q ) & is.null( beta ) )
    {
        if( length( ZI ) == 1 ) ZI = rep( ZI, p )
        if( length( ZI ) != p ) stop( "'ZI', as a vector, must be of length equal to the number of variables, 'ncol( data )'" )
    
        if( is.null( x ) )
        {
            q    = NULL
            beta = NULL
            pii  = NULL
            sample_marginals = vector( "list", p )
            
            cat( paste( c( " MCMC sampling of DW regression parameters ... in progress: \n" ), collapse = "" ) ) 
            
            for( j in 1 : p )
            {
                cat( paste( c(" Marginal regression for node ", j, " " ), collapse = "" ) , "\r" )
                
                xy_j = data.frame( y = data[ , j ] ) 
                
                est_qbeta = BDgraph::bdw.reg( data = xy_j, formula = formula, ZI = ZI[ j ], iter = iter_bdw )
                q[    j ] = est_qbeta $ q.est
                beta[ j ] = est_qbeta $ beta.est
                pii[  j ] = est_qbeta $ pi.est
                sample_marginals[[ j ]] = est_qbeta $ sample
            }
        }
        
        if( !is.null( x ) )
        {
            q    = data * 0
            beta = data * 0
            pii  = vector( length = p )
            sample_marginals = vector( "list", p )
            
            cat( paste( c( " MCMC sampling of DW regression parameters... in progress: \n" ), collapse = "" ) ) 
            
            for( j in 1 : p )
            {
                cat( paste( c(" Marginal regression for node ", j," " ), collapse = "" ) , "\r" )
        
                xy_j = data.frame( x, y = data[ , j ] ) 
                
                est_qbeta   = BDgraph::bdw.reg( data = xy_j, formula = formula, ZI = ZI[ j ], iter = iter_bdw)
                q[    , j ] = est_qbeta $ q.est
                beta[ , j ] = est_qbeta $ beta.est
                pii[    j ] = est_qbeta $ pi.est
                sample_marginals[[ j ]] = est_qbeta $ sample
            }
        }
    }
    
    b      = df.prior
    b_star = b + n
    D      = diag( p )
    Ds     = D + S
    Ts     = chol( solve( Ds ) )
    Ti     = chol( solve( D  ) )   # only for double Metropolis-Hastings algorithms 
    
    g_prior = BDgraph::get_g_prior( g.prior = g.prior, p = p )
    G       = BDgraph::get_g_start( g.start = g.start, g_prior = g_prior, p = p )
    K       = BDgraph::get_K_start( G = G, g.start = g.start, Ts = Ts, b_star = b_star, threshold = threshold )
    
    if( save == TRUE )
    {
        qp1           = ( p * ( p - 1 ) / 2 ) + 1
        string_g      = paste( c( rep( 0, qp1 ) ), collapse = '' )
        sample_graphs = c( rep ( string_g, iter - burnin ) )  # vector of numbers like "10100" 
        graph_weights = c( rep ( 0, iter - burnin ) )         # waiting time for every state
        all_graphs    = c( rep ( 0, iter - burnin ) )         # vector of numbers like "10100"
        all_weights   = c( rep ( 1, iter - burnin ) )         # waiting time for every state		
        size_sample_g = 0
    }else{
        p_links = matrix( 0, p, p )
    }
    
    if( ( verbose == TRUE ) && ( save == TRUE ) && ( p > 50 & iter > 20000 ) )
    {
        cat( "  WARNING: Memory needed to run this function is around " )
        print( ( iter - burnin ) * utils::object.size( string_g ), units = "auto" ) 
    } 
    
    K_hat      = matrix( 0, p, p )
    last_graph = K_hat
    last_K     = K_hat
    
    if( ( is.null( jump ) ) && ( p > 10 & iter > ( 5000 / p ) ) )
        jump = floor( p / 10 )
    
    if( is.null( jump ) ) jump = 1
    
    if( ( p < 10 ) && ( jump > 1 ) )      cat( " WARNING: the value of jump should be 1 " )
    if( jump > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: the value of jump should be smaller " )
    
    if( verbose == TRUE ) 
        cat( paste( c( iter, " MCMC sampling ... in progress: \n" ), collapse = "" ) ) 
    
    bounds = BDgraph::get_bounds_dw( data = data, q = q, beta = beta, pii = pii, n = n, p = p )
    lower_bounds = bounds $ lower_bounds    
    upper_bounds = bounds $ upper_bounds    
    
    # - -  main BDMCMC algorithms implemented in C++ - - - - - - - - - - - - - |
    if( save == TRUE )
    {
        if( ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
        {
            result = .C( "gcgm_dw_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), 
                         K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(data), as.double(lower_bounds), as.double(upper_bounds), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
        {
            counter_all_g  = 0
            
            result = .C( "gcgm_dw_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(data), as.double(lower_bounds), as.double(upper_bounds), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
    }else{
        
        if( ( algorithm == "bdmcmc" ) && ( jump == 1 )  )
        {
            result = .C( "gcgm_dw_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(data), as.double(lower_bounds), as.double(upper_bounds), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
        {
            result = .C( "gcgm_dw_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(data), as.double(lower_bounds), as.double(upper_bounds), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
    }
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    
    K_hat      = matrix( result $ K_hat, p, p, dimnames = list( colnames_data, colnames_data ) ) 
    last_graph = matrix( result $ G    , p, p, dimnames = list( colnames_data, colnames_data ) )
    last_K     = matrix( result $ K    , p, p )
    
    if( save == TRUE )
    {
        if( algorithm == "rjmcmc" ) K_hat = K_hat / ( iter - burnin )		
        size_sample_g = result $ size_sample_g
        sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
        graph_weights = result $ graph_weights[ 1 : size_sample_g ]
        all_graphs    = result $ all_graphs + 1
        all_weights   = result $ all_weights
        if( ( algorithm != "rjmcmc" ) & ( jump != 1 ) )
        { 
            all_weights = all_weights[ 1 : ( result $ counter_all_g ) ]
            all_graphs  = all_graphs[  1 : ( result $ counter_all_g ) ] 
        }
        
        output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = K_hat, 
                       all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph, last_K = last_K, 
                       q.est = q, beta.est = beta, pi.est = pii,
                       data = data, method = "dw" )

    }else{
        p_links = matrix( result $ p_links, p, p, dimnames = list( colnames_data, colnames_data ) ) 
        
        if( ( algorithm == "rjmcmc" ) | ( algorithm == "rj-dmh" ) )
        {
            p_links = p_links / ( iter - burnin )
            K_hat   = K_hat   / ( iter - burnin )
        }
        
        p_links[ lower.tri( p_links ) ] = 0
        
        output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K, 
                       q.est = q, beta.est = beta, pi.est = pii,
                       data = data, method = "dw" )
    }
    
    if( !is.null( sample_marginals ) ) 
        output $ sample_marginals = sample_marginals
    
    class( output ) = "bdgraph"
    return( output )   
}
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
