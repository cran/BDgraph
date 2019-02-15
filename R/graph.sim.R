## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
#     Copyright (C) 2012 - 2019  Reza Mohammadi                                                    |
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
    if( ( sum( prob < 0 ) + sum( prob > 1 ) ) != 0 ) stop( "'prob' must be between 0 and 1" )
 
    G <- matrix( 0, p, p )
    
    # - - build the graph structure - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    if( graph == "random" )
    {
        if( is.null( size ) )
        {
            G[ upper.tri( G ) ] <- stats::rbinom( p * ( p - 1 ) / 2, 1, prob )
        }else{
            if( size < 0 | size > p * ( p - 1 ) / 2 )  stop( "Graph size must be between zero and p*(p-1)/2" )
            
            smp <- sample( 1 : ( p * ( p - 1 ) / 2 ), size, replace = FALSE )
            G[ upper.tri( G ) ][smp] <- 1
        }
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
            #if( !is.null( size ) ) class <- length( size ) else class <- max( 2, ceiling( p / 20 ) )
        }
        
        g.large <- p %% class
        g.small <- class - g.large
        n.small <- floor( p / class )
        n.large <- n.small + 1
        vp      <- c( rep( n.small, g.small ), rep( n.large, g.large ) )
        
        if( is.null( size ) )
        {
            if( length( prob ) != class ) prob = rep( prob, class )
            
            for( i in 1 : class )
            {
                tmp <- if( i == 1 ) ( 1 : vp[ 1 ] ) else ( ( sum( vp[ 1 : ( i - 1 ) ] ) + 1 ) : sum( vp[ 1 : i ] ) )
                gg                <- matrix( 0, vp[ i ], vp[ i ] )
                gg[ upper.tri( gg ) ] <- stats::rbinom( vp[ i ] * ( vp[ i ] - 1 ) / 2, 1, prob[ i ] )
                G[ tmp, tmp ]       <- gg
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
    }
    
    if( graph == "scale-free" )
    {
        resultGraph = .C( "scale_free", G = as.integer(G), as.integer(p), PACKAGE = "BDgraph" )
        G = matrix( resultGraph $ G, p, p ) 
        
        #j = sample( 1:p, 1 )
        #for( i in ( c( 1:p )[ -j ] ) ) { G[ i, j ] = 1; G[ j, i ] = 1 }
    }
    
    if( ( graph == "lattice" ) | ( graph == "grid" ) )
    {
        if( is.null( size ) )
        {
            length_row = round( sqrt( p ) )
            length_col = round( sqrt( p ) )
        }else{
            if( length( size ) == 1 )
            {
                length_row = size
                length_col = size
            }else{
                length_row = size[ 1 ]
                length_col = size[ 2 ]
            } 
        }
        
        for( row in 1:length_row )
        {
            for( col in 1:length_col )
            {
                if( ( row != length_row ) & ( col != length_col ) )
                    G[ col + ( row - 1 ) * length_col, c( col + ( row - 1 ) * length_col + 1, col + row * length_col ) ] = 1
                if( ( row == length_row ) & ( col != length_col ) )
                    G[ col + ( row - 1 ) * length_col, col + ( row - 1 ) * length_col + 1 ] = 1
                if( ( row != length_row ) & ( col == length_col ) )
                    G[ col + ( row - 1 ) * length_col, col + row * length_col ] = 1
            }
        }
    }
    
    if( graph == "hub" )
    {
        if( is.null( size ) ) size = ceiling( p / 20 ) 

        hub = sample( 1:p, size = size, replace = FALSE )
        
        for( i in 1:size )
        {
            G[ hub[ i ],  ] <- 1
            G[ , hub[ i ] ] <- 1
        }
    }

    if( graph == "star" )
    {
        hub = sample( 1:p, size = 1, replace = FALSE )
        G[ hub,  ] <- 1
        G[ , hub ] <- 1
    }
    
    if( graph == "circle" )
    {
        if( p < 2 ) stop( "For 'circle' graph, 'p' must be more than 2" )
        
        G         <- stats::toeplitz( c( 0, 1, rep( 0, p - 2 ) ) )
        G[ 1, p ] <- 1
    }
    
    G[ lower.tri( G, diag = T ) ] = 0
    G = G + t( G )

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

