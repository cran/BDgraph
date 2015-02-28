# To select the graph in which the edge posterior probabilities are more than "cut" value
# OR if cut is NULL to select the best graph (graph with the highest posterior probability) 
select = function ( output, cut = NULL, vis = FALSE )
{
	if ( is.null(cut) )
	{
		sampleGraphs <- output $ sampleGraphs
		graphWeights <- output $ graphWeights
		p            <- nrow( output $ lastGraph )
		prob.G       <- graphWeights / sum( graphWeights )
		max.prob.G   <- which( prob.G == max( prob.G ) )
		
		if( length(max.prob.G) > 1 ) max.prob.G <- max.prob.G[1] 
		
		gi        <- sampleGraphs[max.prob.G]
		gv        <- c( rep( 0, p * ( p - 1 ) / 2 ) )
		gv[ which( unlist( strsplit( as.character(gi), "" ) ) == 1 ) ] <- 1

		dimlab   <- colnames( output $ lastGraph )
		if ( is.null(dimlab) ) dimlab <- as.character(1 : p)
		
		graphi <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	

		graphi[ upper.tri(graphi) ] <- gv
	} 
	else 
	{
		if ( ( cut < 0 ) || ( cut > 1 ) )   stop( "Value of 'cut' should be between zero and one." )
		prob  <- as.matrix( phat( output ) )
		prob[ prob > cut ] = 1
		prob[ prob <= cut ]  = 0
		graphi = prob
	}
	
	if ( vis )
	{
		G <- graph.adjacency( graphi, mode = "undirected", diag = FALSE )
		if ( p < 20 ) sizev = 15 else sizev = 2

		if ( is.null(cut) )
		{
			plot.igraph( G, layout = layout.circle, main = "Graph with highest posterior probability", sub = paste( c( "Posterior probability = ", round( max( prob.G ), 4) ), collapse = "" ),
			            vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black'  )
		} 
		else 
		{
			plot.igraph( G, layout = layout.circle, main = paste( c( "Graph with links posterior probabilities = ",  cut ), collapse = "" ), vertex.color = "white", vertex.size = sizev, vertex.label.color = 'black' )
		}
	}

	return( Matrix( graphi + t(graphi), sparse = TRUE ) )
}
     
