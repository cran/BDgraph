# computing probability of all links of the graph
phat = function( x, round = 3 )
{
	phat = x $ phat

	if( is.null( phat ) )
	{
		sampleGraphs <- x $ sampleGraphs
		graphWeights <- x $ graphWeights
		p            <- nrow( x $ lastGraph )
		pvec         <- c( rep( 0, p * ( p - 1 ) / 2) )
	   
		for ( i in 1 : length( sampleGraphs ) )
		{
			inp       <- which( unlist( strsplit( as.character( sampleGraphs[i] ), "" ) ) == 1 )
			pvec[inp] <- pvec[inp] + graphWeights[i]
		}
		
		dimlab <- colnames( x $ lastGraph ) # lastG
		phat   <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )
		phat[ upper.tri(phat) ] <- pvec / sum( graphWeights )
	}
	
	return( Matrix( round( phat, round ) ) )
}
      
