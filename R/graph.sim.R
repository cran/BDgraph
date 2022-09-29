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
#     Graph generator                                                          |
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

graph.sim = function( p = 10, graph = "random", prob = 0.2, size = NULL, 
                      class = NULL, vis = FALSE, rewire = 0.05 )
{
    if( p < 2 ) stop( "'p' must be more than 1" )
    if( ( prob < 0 ) || ( prob > 1 )    ) stop( "'prob' must be between ( 0, 1 )" )
    if( ( rewire < 0 ) | ( rewire > 1 ) ) stop( "Value of 'rewire' must be between ( 0, 1 )" )

    G <- matrix( 0, p, p )
    
    # - - build the graph structure - - - - - - - - - - - - - - - - - - - - -  |
    if( ( graph == "random" ) | ( graph == "Random" ) )
    {
        if( is.null( size ) )
        {
            G[ upper.tri( G ) ] <- stats::rbinom( p * ( p - 1 ) / 2, 1, prob )
        }else{
            if( ( size < 0 ) | ( size > p * ( p - 1 ) / 2 ) )  stop( "'size' must be between ( 0,  p * ( p - 1 ) / 2 )" )
            
            smp <- sample( 1 : ( p * ( p - 1 ) / 2 ), size, replace = FALSE )
            G[ upper.tri( G ) ][ smp ] <- 1
        }
    }

    if( ( graph == "scale-free" ) | ( graph == "Scale-free" ) )
    {
        resultGraph = .C( "scale_free", G = as.integer( G ), as.integer( p ), PACKAGE = "BDgraph" )
        G = matrix( resultGraph $ G, p, p ) 
        
        #j = sample( 1:p, 1 )
        #for( i in ( c( 1:p )[ -j ] ) ) { G[ i, j ] = 1; G[ j, i ] = 1 }
    }
    
    if( ( graph == "cluster" ) | ( graph == "Cluster" ) )
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
            if( class != length( size ) )  stop( "Number of graph sizes is not match with number of clusters" )
            if( ( sum( size ) < 0 ) | ( sum( size ) > p * ( p - 1 ) / 2 ) ) stop( "Total graph sizes must be between ( 0, p * ( p - 1 ) / 2 )" )
            
            for( i in 1 : class )
            {
                tmp <- if( i == 1 ) ( 1 : vp[ 1 ] ) else ( ( sum( vp[ 1 : ( i - 1 ) ] ) + 1 ) : sum( vp[ 1 : i ] ) )
                gg  <- matrix( 0, vp[ i ], vp[ i ] )
                smp <- sample( 1 : ( vp[ i ] * ( vp[ i ] - 1 ) / 2 ), size[ i ], replace = FALSE )
                gg[ upper.tri( gg ) ][ smp ] <- 1
                G[ tmp, tmp ]            <- gg
            }
        }
    }
    
    if( ( graph == "hub" ) | ( graph == "Hub" ) )
    {
        if( is.null( size ) ) size = ceiling( p / 20 ) 
        if( ( size < 0 ) | ( size > ( p - 1 ) ) )  stop( "'size' must be between ( 0, p - 1 ), for option 'graph = \"hub\"'" )

        hub = sample( 1:p, size = size, replace = FALSE )
        
        for( i in 1:size )
        {
            G[ hub[ i ],  ] <- 1
            G[ , hub[ i ] ] <- 1
        }
    }

    if( ( graph == "star" ) | ( graph == "Star" ) )
    {
        hub = sample( 1:p, size = 1, replace = FALSE )
        G[ hub,     ] <- 1
        G[    , hub ] <- 1
    }
    
    if( ( graph == "circle" ) | ( graph == "Circle" ) )
    {
        if( p < 3 ) stop( "'p' must be more than 2, for option 'graph = \"circle\"'" )
        
        G         <- stats::toeplitz( c( 0, 1, rep( 0, p - 2 ) ) )
        G[ 1, p ] <- 1
    }

    if( ( graph == "smallworld" ) | ( graph == "Smallworld" ) | ( graph == "small-world" ) | ( graph == "Small-world" ) )
    {
        G_igraph = igraph::sample_smallworld( dim = 1, # One dimension
                                              size = p, # Number of variables
                                              nei = round( size / p ), # Neighborhood
                                              p = rewire )   
        
        G = as.matrix( igraph::as_adj( G_igraph ) ) # Rewiring probability
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
    
    G[ lower.tri( G, diag = TRUE ) ] = 0
    G = G + t( G )

    # - - graph visualization - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|
    if( vis == TRUE )
        BDgraph::plot.graph( G, main = "Graph structure" )
    
    class( G ) <- "graph"
    return( G )
}
   
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
# plot for class "graph" from graph.sim function
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

plot.graph = function( x, cut = 0.5, 
                       mode = "undirected", diag = FALSE, main = NULL, 
                       layout = igraph::layout_with_fr, 
                       vertex.size = 2, 
                       vertex.color = "orange", 
                       vertex.frame.color = "orange", 
                       vertex.label = NULL,
                       vertex.label.dist = 0.5, 
                       vertex.label.color = "blue", 
                       edge.color = "lightblue", ... )
{
    graph = BDgraph::get_graph( x, cut = cut )
    
    if( is.null( vertex.label ) ) vertex.label = colnames( graph )
    
    graph_ig <- igraph::graph.adjacency( graph, mode = mode, diag = diag )
    
    igraph::plot.igraph( graph_ig, 
                         main = main, 
                         layout = layout, 
                         vertex.size = vertex.size,
                         vertex.color = vertex.color, 
                         vertex.frame.color = vertex.frame.color,
                         vertex.label = vertex.label,
                         vertex.label.dist = vertex.label.dist, 
                         vertex.label.color = vertex.label.color,
                         edge.color = edge.color, ... )
}		
    
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
