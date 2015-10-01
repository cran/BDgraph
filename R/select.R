# To select the graph in which the edge posterior probabilities are more than "cut" value
# OR if cut is NULL to select the best graph (graph with the highest posterior probability) 
select = function( x, cut = 0.5, vis = FALSE )
{
	phat = x $ phat
	p    = nrow( x $ lastGraph )
  
	if( is.null( phat ) )
	{
		if( is.null( cut ) )
		{
			sampleGraphs <- x $ sampleGraphs
			graphWeights <- x $ graphWeights
			
			indG_max     <- sampleGraphs[ which( graphWeights == max( graphWeights ) )[1] ]
			vec_G        <- c( rep( 0, p * ( p - 1 ) / 2 ) )
			vec_G[ which( unlist( strsplit( as.character( indG_max ), "" ) ) == 1 ) ] <- 1

			dimlab       <- colnames( x $ lastGraph )
			selected_G   <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	
			selected_G[ upper.tri(selected_G) ] <- vec_G
		} 
		else 
		{
			if ( ( cut < 0 ) || ( cut > 1 ) ) stop( "Value of 'cut' should be between zero and one." )
			phat                = as.matrix( phat( x ) )
			phat[ phat > cut ]  = 1
			phat[ phat <= cut ] = 0
			selected_G          = phat
		}
	}
	else
	{
		if( ( cut < 0 ) || ( cut > 1 ) ) stop( "Value of 'cut' should be between zero and one." )
		selected_G                = 0 * phat
		selected_G[ phat > cut ]  = 1
		selected_G[ phat <= cut ] = 0
	}
		
	if( vis )
	{
		G <- graph.adjacency( selected_G, mode = "undirected", diag = FALSE )
		if( p < 20 ) sizev = 15 else sizev = 2

		if( is.null(cut) )
		{
			plot.igraph( G, layout = layout.circle, main = "Graph with highest posterior probability", sub = paste( c( "Posterior probability = ", round( max( graphWeights ) / sum( graphWeights ), 4) ), collapse = "" ),
			            vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black'  )
		} 
		else 
		{
			plot.igraph( G, layout = layout.circle, main = paste( c( "Graph with links posterior probabilities > ",  cut ), collapse = "" ), vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		}
	}

	return( Matrix( selected_G, sparse = TRUE ) )
}
       
