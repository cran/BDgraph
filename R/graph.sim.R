## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2018  Reza Mohammadi                                                    |
#                                                                                                  |
#     This file is part of BDgraph package.                                                        |
#                                                                                                  |
#     BDgraph is free software: you can redistribute it and/or modify it under                     |
#     the terms of the GNU General Public License as published by the Free                         |
#     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                    |
#                                                                                                  |
#     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Graph generator                                                                              |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

graph.sim = function( p = 10, graph = "random", prob = 0.2, size = NULL, class = NULL, vis = FALSE )
{
    if( p < 2 ) stop( "'p' must be more than 1" )
    if( prob < 0 | prob > 1 ) stop( "'prob' must be between zero and one" )

    # - - build the graph structure - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    if( graph == "random" )
    {
        G <- matrix( 0, p, p )
        
        if( is.null( size ) )
        {
            G[ upper.tri( G ) ] <- stats::rbinom( p * ( p - 1 ) / 2, 1, prob )
        }else{
            if( size < 0 | size > p * ( p - 1 ) / 2 )  stop( "Graph size must be between zero and p*(p-1)/2" )
            
            smp <- sample( 1 : ( p * ( p - 1 ) / 2 ), size, replace = FALSE )
            G[ upper.tri( G ) ][smp] <- 1
        }
        
        G <- G + t( G )
    }
    
    if( graph == "cluster" )
    {
        # partition variables
        if( is.null( class ) )
        { 
            #class = NULL
            if( !is.null( size ) )   class = length( size )
            if( length( prob ) > 1 ) class = length( prob )
            if( is.null( class ) )   class = max( 2, ceiling( p / 20 ) )
            
            if( !is.null( size ) )
            {
                class <- length( size )
            }else{
                class <- max( 2, ceiling( p / 20 ) )
            }
        }
        
        g.large <- p %% class
        g.small <- class - g.large
        n.small <- floor( p / class )
        n.large <- n.small + 1
        vp      <- c( rep( n.small, g.small ), rep( n.large, g.large ) )
        
        G       <- matrix( 0, p, p )
        
        if( is.null( size ) )
        {
            # if( prob < 0 | prob > 1 ) stop( " 'prob' must be between zero and one" )
            
            for( i in 1 : class )
            {
                tmp <- if( i == 1 ) ( 1 : vp[1] ) else ( ( sum( vp[1 : (i-1)] ) + 1 ) : sum( vp[1:i] ) )
                gg                <- matrix( 0, vp[i], vp[i] )
                gg[upper.tri(gg)] <- stats::rbinom( vp[i] * ( vp[i] - 1 ) / 2, 1, prob )
                G[tmp, tmp]       <- gg
            }
        }else{
            if( class != length( size ) )  stop( " Number of graph sizes is not match with number of clusters" )
            if( sum( size ) < 0 | sum( size ) > p * ( p - 1 ) / 2 ) stop( " Total graph sizes must be between zero and p*(p-1)/2" )
            
            for( i in 1 : class )
            {
                tmp <- if( i == 1 ) ( 1 : vp[1] ) else ( ( sum( vp[1 : (i-1)] ) + 1 ) : sum( vp[1:i] ) )
                gg  <- matrix( 0, vp[i], vp[i] )
                smp <- sample( 1 : ( vp[i] * (vp[i] - 1) / 2 ), size[i], replace = FALSE )
                gg[upper.tri(gg)][smp] <- 1
                G[tmp, tmp]            <- gg
            }
        }
        
        G <- G + t( G )	   
    }
    
    if( graph == "hub" )
    {
        # partition variables
        if( is.null( class ) ) class <- ceiling( p / 20 ) 
        
        # partition variables into groups
        g.large <- p %% class
        g.small <- class - g.large
        n.small <- floor( p / class )
        n.large <- n.small + 1
        g.list  <- c( rep( n.small, g.small ), rep( n.large, g.large ) )
        g.ind   <- rep( c( 1:class ), g.list )
        
        G <- matrix( 0, p, p )
        
        for( i in 1:class )
        {
            tmp              <- which( g.ind == i )
            G[ tmp[1], tmp ] <- 1
            G[ tmp, tmp[1] ] <- 1
        }
    }
    
    if( graph == "circle" )
    {
        if( p < 2 ) stop( "For 'circle' graph, 'p' must be more than 2" )
        
        G         <- stats::toeplitz( c( 0, 1, rep( 0, p - 2 ) ) )
        G[ 1, p ] <- 1
        G[ p, 1 ] <- 1
    }
    
    if( graph == "scale-free" )
    {
        G = matrix( 0, p, p )
        resultGraph = .C( "scale_free", G = as.integer(G), as.integer(p), PACKAGE = "BDgraph" )
        G = matrix( resultGraph $ G, p, p ) 
        
        j = sample( 1:p, 1 )
        
        for( i in ( c( 1:p )[ -j ] ) )
        {
            G[ i, j ] = 1
            G[ j, i ] = 1
        }
    }
    
    if( ( graph == "lattice" ) | ( graph == "grid" ) )
    {
        G = matrix( 0, p, p )
        sp = floor( sqrt( p ) )
        
        for( row in 1:( sp - 1 ) )
        {
            for( col in 1:( sp - 1 ) )
            {
                node = ( row - 1 ) * sp + col
                G[ node, node + 1  ] = 1
                G[ node, node + sp ] = 1
            }
            
            node = sp * row
            G[ node, node + sp ] = 1
            
            node = sp * ( sp - 1 ) + row
            G[ node, node + 1 ] = 1
        }        
    }
        
    # - - graph visualization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    if( vis == TRUE )
    {
        graph_ig <- igraph::graph.adjacency( G, mode = "undirected", diag = FALSE )
        
        if( p < 20 ) size = 10 else size = 2
        igraph::plot.igraph( graph_ig, layout = igraph::layout.circle, main = "Graph structure", 
                     vertex.color = "white", vertex.size = size, vertex.label.color = 'black' )
    }
    
    class( G ) <- "graph"
    return( G )
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# plot for class "graph" from graph.sim function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
plot.graph = function( x, main = NULL, layout = layout.circle, ... )
{
    true_graph = as.matrix( x )
    if( is.null( main ) ) main = "Graph structure"
    g_igraph <- igraph::graph.adjacency( true_graph, mode = "undirected", diag = FALSE )
    
    igraph::plot.igraph( g_igraph, main = main, layout = layout, ... )
}		

