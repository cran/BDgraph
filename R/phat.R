# computing probability of all links of the graph
phat = function( output, round = 3 )
{
	sampleGraphs <- output $ sampleGraphs
	graphWeights <- output $ graphWeights
	p            <- nrow( output $ lastGraph )
	pvec         <- c( rep( 0, p * ( p - 1 ) / 2) )
   
	for ( i in 1 : length( sampleGraphs ) )
	{
		inp       <- which( unlist( strsplit( as.character( sampleGraphs[i] ), "" ) ) == 1 )
		pvec[inp] <- pvec[inp] + graphWeights[i]
	}
	
	dimlab <- colnames( output $ lastGraph ) # lastG
	if ( is.null( dimlab ) ) dimlab <- as.character( 1 : p )
	
	phat   <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )
	phat[ upper.tri(phat) ] <- pvec / sum( graphWeights )

	return( Matrix( round( phat, round ) ) )
}
    
