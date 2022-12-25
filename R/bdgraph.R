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
#     Main function of BDgraph package: BDMCMC algorithm for graphical models  |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

bdgraph = function( data, n = NULL, method = "ggm", algorithm = "bdmcmc", iter = 5000, 
                    burnin = iter / 2, not.cont = NULL, g.prior = 0.2, df.prior = 3, 
                    g.start = "empty", jump = NULL, save = FALSE, 
                    cores = NULL, threshold = 1e-8, verbose = TRUE, nu = 1 )
{
    if( df.prior < 3  ) stop( "'prior.df' must be >= 3" )
    if( iter < burnin ) stop( "'iter' must be higher than 'burnin'" )
    
    burnin <- floor( burnin )
    
    if( is.numeric( verbose ) )
    {
        if( ( verbose < 1 ) | ( verbose > 100 ) ) 
            stop( "'verbose' (for numeric case) must be between ( 1, 100 )" )
        
        trace_mcmc = floor( verbose )
        verbose = TRUE
        
    }else{
        trace_mcmc = ifelse( verbose == TRUE, 10, iter + 1000 )
    }

    list_S_n_p = BDgraph::get_S_n_p( data = data, method = method, n = n, not.cont = not.cont )
    
    S      = list_S_n_p $ S
    n      = list_S_n_p $ n
    p      = list_S_n_p $ p
    method = list_S_n_p $ method
    data   = list_S_n_p $ data
    colnames_data = list_S_n_p $ colnames_data
    
    if( ( is.null( cores ) ) & ( p < 16 ) ) 
        cours = 1
        
    cores = BDgraph::get_cores( cores = cores, verbose = verbose )
    
    if( method == "gcgm" )
    {
        not.cont = list_S_n_p $ not.cont
        R        = list_S_n_p $ R
        Z        = list_S_n_p $ Z
        gcgm_NA  = list_S_n_p $ gcgm_NA
    }
 
    if( method == "tgm" )
    {
        tu = stats::rgamma( n, shape = nu / 2, rate = nu / 2 )
        
        mu = tu %*% data / sum( tu )
        
        S = matrix( 0, p, p )
        
        for( i in 1:n )	
        { 
            d_mu = data[ i, , drop = FALSE ] - mu
            S    = S + tu[ i ] * t( d_mu ) %*% d_mu	
        }
    }
       
    b      = df.prior
    b_star = b + n
    D      = diag( p )
    Ds     = D + S
    Ts     = chol( solve( Ds ) )
    Ti     = chol( solve( D ) )   # only for double Metropolis-Hastings algorithms 
    
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
        cat( "  WARNING: Memory needs to run this function is around: " )
        print( ( iter - burnin ) * utils::object.size( string_g ), units = "auto" ) 
    } 
    
    K_hat      = matrix( 0, p, p )
    last_graph = K_hat
    last_K     = K_hat
    
    if( ( is.null( jump ) ) && ( p > 10 & iter > ( 5000 / p ) ) )
        jump = floor( p / 10 )
    
    if( is.null( jump ) ) jump = 1
    
    if( ( p < 10 ) && ( jump > 1 ) )      cat( " WARNING: the value of jump should be 1. " )
    if( jump > min( p, sqrt( p * 11 ) ) ) cat( " WARNING: the value of jump should be smaller. " )
    
    if( verbose == TRUE ) 
        cat( paste( c( iter, " MCMC sampling ... in progress: \n" ), collapse = "" ) ) 
    
    # - -  main BDMCMC algorithms implemented in C++ - - - - - - - - - - - - - |
    if( save == TRUE )
    {
        if( ( method == "tgm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
        {
            result = .C( "tgm_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), 
                         as.double(data), as.integer(n), as.double(nu), as.double(mu), as.double(tu),
                         PACKAGE = "BDgraph" )
        }

        if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
        {
            result = .C( "ggm_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold), 
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
        {
            counter_all_g = 0
            result = .C( "ggm_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold), 
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
        {
            result = .C( "ggm_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold), 
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
        {
            not_continuous = not.cont
            counter_all_g  = 0
            
            result = .C( "gcgm_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "rjmcmc" ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }	
        
        # for Double Metropolis-Hasting 
        if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( jump == 1 ) )
        {
            result = .C( "ggm_DMH_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold), 
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( jump != 1 ) )
        {
            counter_all_g = 0
            result = .C( "ggm_DMH_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold), 
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "ggm" ) && ( algorithm == "rj-dmh" ) )
        {
            result = .C( "ggm_DMH_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold), 
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( jump == 1 ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_DMH_bdmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( jump != 1 ) )
        {
            not_continuous   = not.cont
            counter_all_g = 0
            
            result = .C( "gcgm_DMH_bdmcmc_map_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g), counter_all_g = as.integer(counter_all_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "rj-dmh" ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_DMH_rjmcmc_map", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         all_graphs = as.integer(all_graphs), all_weights = as.double(all_weights), K_hat = as.double(K_hat), 
                         sample_graphs = as.character(sample_graphs), graph_weights = as.double(graph_weights), size_sample_g = as.integer(size_sample_g),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }	
        
    }else{

        if( ( method == "tgm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 )  )
        {
            
            # double D[], double data[], int *n, double *nu, double mu[], double tu[]
             result = .C( "tgm_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.integer(trace_mcmc), 
                         as.double(D), as.double(data), as.integer(n), as.double(nu), as.double(mu), as.double(tu),
                         PACKAGE = "BDgraph" )
        }
                
        if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 )  )
        {
            result = .C( "ggm_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
        {
            result = .C( "ggm_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }		
        
        if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
        {
            result = .C( "ggm_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.integer(p_links),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( jump == 1 )  )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) && ( jump != 1 ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "rjmcmc" ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.integer(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }	
        
        # for Double Metropolis-Hasting 
        if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( jump == 1 )  )
        {
            result = .C( "ggm_DMH_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "ggm" ) && ( algorithm == "bd-dmh" ) && ( jump != 1 ) )
        {
            result = .C( "ggm_DMH_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }		
        
        if( ( method == "ggm" ) && ( algorithm == "rj-dmh" ) )
        {
            result = .C( "ggm_DMH_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold), 
                         K_hat = as.double(K_hat), p_links = as.integer(p_links),
                         as.integer(b), as.integer(b_star), as.double(Ds), as.double(D), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( jump == 1 )  )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_DMH_bdmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "bd-dmh" ) && ( jump != 1 ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_DMH_bdmcmc_ma_multi_update", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.double(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(jump), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }
        
        if( ( method == "gcgm" ) && ( algorithm == "rj-dmh" ) )
        {
            not_continuous = not.cont
            
            result = .C( "gcgm_DMH_rjmcmc_ma", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(g_prior), as.double(Ts), as.double(Ti), K = as.double(K), as.integer(p), as.double(threshold),
                         as.double(Z), as.integer(R), as.integer(not_continuous), as.integer(n), as.integer(gcgm_NA),
                         K_hat = as.double(K_hat), p_links = as.integer(p_links),
                         as.integer(b), as.integer(b_star), as.double(D), as.double(Ds), as.integer(trace_mcmc), PACKAGE = "BDgraph" )
        }	
    }
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|

    K_hat      = matrix( result $ K_hat, p, p, dimnames = list( colnames_data, colnames_data ) ) 
    last_graph = matrix( result $ G    , p, p, dimnames = list( colnames_data, colnames_data ) )
    last_K     = matrix( result $ K    , p, p )
    
    if( save == TRUE )
    {
        if( ( algorithm == "rjmcmc" ) | ( algorithm == "rj-dmh" ) )
            K_hat = K_hat / ( iter - burnin )		
        
        size_sample_g = result $ size_sample_g
        sample_graphs = result $ sample_graphs[ 1 : size_sample_g ]
        graph_weights = result $ graph_weights[ 1 : size_sample_g ]
        all_graphs    = result $ all_graphs + 1
        all_weights   = result $ all_weights
        if( ( jump != 1 ) & ( algorithm != "rjmcmc" ) & ( algorithm != "rj-dmh" ) )
        { 
            all_weights = all_weights[ 1 : ( result $ counter_all_g ) ]
            all_graphs  = all_graphs[  1 : ( result $ counter_all_g ) ] 
        }
        
        output = list( sample_graphs = sample_graphs, graph_weights = graph_weights, K_hat = K_hat, 
                       all_graphs = all_graphs, all_weights = all_weights, last_graph = last_graph, last_K = last_K,
                       data = data, method = method )
    }else{
        p_links = matrix( result $ p_links, p, p, dimnames = list( colnames_data, colnames_data ) ) 
        
        if( ( algorithm == "rjmcmc" ) | ( algorithm == "rj-dmh" ) )
        {
            p_links = p_links / ( iter - burnin )
            K_hat   = K_hat   / ( iter - burnin )
        }
        
        p_links[ lower.tri( p_links ) ] = 0
        output = list( p_links = p_links, K_hat = K_hat, last_graph = last_graph, last_K = last_K,
                       data = data, method = method )
    }
    
    class( output ) = "bdgraph"
    return( output )   
}
  
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Summary for "bdgraph" object                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
summary.bdgraph = function( object, round = 2, vis = TRUE, ... )
{
    K_hat      = object $ K_hat
    p_links    = object $ p_links
	if( is.null( p_links ) ) p_links = BDgraph::plinks( object )

    selected_g = BDgraph::select( p_links, cut = 0.5 )    

	if( vis == TRUE )
	{
		if( !is.null( object $ graph_weights ) ) 
			op = graphics::par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.1, 0.1, 0.1, 0.1 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 
		
		# - - - plot selected graph
		sub_g = "Graph with edge posteriors > 0.5"
		BDgraph::plot.graph( selected_g, main = "Selected Graph", sub = sub_g, ... )

		if( !is.null( object $ graph_weights ) )
		{
		    sample_graphs = object $ sample_graphs
		    graph_weights = object $ graph_weights
		    sum_gWeights  = sum( graph_weights )
		    
		    # - - - plot posterior distribution of graph
		    graph_prob = graph_weights / sum_gWeights
			graphics::plot( x = 1 : length( graph_weights ), y = graph_prob, type = "h", col = "gray60", 
			                main = "Graph Posteriors", ylab = "Pr(graph|data)", 
			                xlab = "Graph", ylim = c( 0, max( graph_prob ) ) )
			
			# - - - plot posterior distribution of graph size
			sizesample_graphs = sapply( sample_graphs, function( x ) length( which( unlist( strsplit( as.character( x ), "" ) ) == 1 ) ) )
			xx       <- unique( sizesample_graphs )
			weightsg <- vector()

			for( i in 1 : length( xx ) ) 
			    weightsg[ i ] <- sum( graph_weights[ which( sizesample_graphs == xx[ i ] ) ] )

			prob_zg = weightsg / sum_gWeights
			graphics::plot( x = xx, y = prob_zg, type = "h", col = "gray10",
			                main = "Graphs Size's Posteriors", 
			                ylab = "Pr(graph size|data)", xlab = "Graph Size",
			                ylim = c( 0, max( prob_zg ) ) )

			# - - - plot trace of graph size
			all_graphs     = object $ all_graphs
			sizeall_graphs = sizesample_graphs[ all_graphs ]
			  
			graphics::plot( x = 1 : length( all_graphs ), sizeall_graphs, type = "l", col = "lightblue", 
			                main = "Trace of Graph Size", ylab = "Graph size", xlab = "Iteration" )

			graphics::par( op )
		}
	}
	
	if( is.null( K_hat ) )			  
		return( list( selected_g = selected_g, p_links = round( p_links, round ) ) )
	else
		return( list( selected_g = selected_g, p_links = round( p_links, round ), K_hat = round( K_hat, round ) ) )
}  
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Plot function for "bdgraph" object                                        |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
plot.bdgraph = function( x, cut = 0.5, number.g = NULL, 
                         main = NULL,
                         layout = igraph::layout_with_fr, 
                         vertex.size = 2, 
                         vertex.color = "orange", 
                         vertex.frame.color = "orange", 
                         vertex.label = NULL,
                         vertex.label.dist = 0.5, 
                         vertex.label.color = "blue", 
                         edge.color = "lightblue", ... )
{
 	if( is.null( number.g ) )
	{
 	    sub = paste0( "Edge posterior probability = ", cut )
 	    
	    BDgraph::plot.graph( x, cut = cut, sub = sub, 
	                         main = main, 
                             layout = layout, 
                             vertex.size = vertex.size,
                             vertex.color = vertex.color, 
                             vertex.frame.color = vertex.frame.color,
	                         vertex.label = vertex.label,
	                         vertex.label.dist = vertex.label.dist, 
                             vertex.label.color = vertex.label.color,
                             edge.color = edge.color, ... )
	}else{
	    
	    if( is.null( x $ all_graphs ) ) stop( "'x' must be an object of 'bdgraph()' function with option 'save = TRUE'" )
	    
	    sample_graphs = x $ sample_graphs
	    graph_weights = x $ graph_weights
	    prob_G        = graph_weights / sum( graph_weights )
	    sort_prob_G   = sort( prob_G, decreasing = TRUE )
	    
	    p             = nrow( x $ last_graph )
	    label         = colnames( x $ last_graph )
	    list_G        = replicate( number.g, matrix( 0, p, p, dimnames = list( label, label ) ), simplify = FALSE )
	    vec_G         = c( rep( 0, p * ( p - 1 ) / 2 ) )
	    
	    if( number.g == 2 ) op <- graphics::par( mfrow = c( 1, 2 ), pty = "s" )
	    if( number.g > 2 & number.g < 7 )  op <- graphics::par( mfrow = c( 2, number.g %% 2 + trunc( number.g / 2 ) ), pty = "s" )
	    
	    for( i in 1 : number.g )
	    {
	        if( number.g > 6 ) grDevices::dev.new()  
	        
	        indG_i <- sample_graphs[ which( prob_G == sort_prob_G[i] )[1] ]
	        vec_G  <- 0 * vec_G
	        vec_G[ which( unlist( strsplit( as.character(indG_i), "" ) ) == 1 ) ] <- 1
	        list_G[[i]][ upper.tri( list_G[[i]] ) ] <- vec_G
	        
	        main = ifelse( i == 1, "Graph with highest probability", paste( c( i, "th graph" ), collapse = "" ) )
	        sub  = paste( c( "Posterior probability = ", round( sort_prob_G[i], 6 ) ), collapse = "" )
	        
	        BDgraph::plot.graph( list_G[[i]], main = main, sub = sub, 
	                             layout = layout, 
                                 vertex.size = vertex.size,
                                 vertex.color = vertex.color, 
                                 vertex.frame.color = vertex.frame.color,
                                 vertex.label.dist = vertex.label.dist, 
                                 vertex.label.color = vertex.label.color,
                                 edge.color = edge.color, ... )
	    }
	    
	    if( number.g > 1 & number.g < 7 ) graphics::par( op )
    }
}
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    Print function for "bdgraph" object                                       |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
print.bdgraph = function( x, ... )
{
	p_links = x $ p_links
	
	if( is.null( p_links ) ) 
	    p_links = BDgraph::plinks( x )
	
	selected_g = BDgraph::select( p_links, cut = 0.5 )

	cat( paste( "\n Adjacency matrix of selected graph \n" ), fill = TRUE )
	print( selected_g )
	
    cat( paste( "\n Edge posterior probability of the links \n" ), fill = TRUE )
    print( round( p_links, 2 ) )
} 
     
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#    predict function for "bdgraph" object                                       |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

predict.bdgraph = function( object, iter = 1, ... )
{
    method = object $ method
    data   = object $ data
    
    n_data = nrow( data )
    p      = ncol( data )

    K = object $ K_hat
    
    if( is.null( K ) )
    {
		if( isSymmetric( data ) )
		{
			S = data
		}else{
 			S = t( data ) %*% data
		}
        
        G = BDgraph::select( bdgraph.obj = object )
        
        sample_K = BDgraph::rgwish( n = 500, adj = G, b = 3 + n_data, D = diag( p ) + S )
        
        K = 0 * G
        for( i in 1:dim( sample_K )[3] )
            K = K + sample_K[[i]]
        
        K = K / dim( sample_K )[3]
    }
    
    sigma = solve( K )
    
    Z = BDgraph::rmvnorm( n = iter, mean = 0, sigma = sigma )
    
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
        
        Z = tmvtnorm::rtmvnorm( n = iter, mean = mean, 
                               sigma = sigma, lower = rep( -5, length = p ), 
                               upper = rep( 5, length = p ) )
        
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














