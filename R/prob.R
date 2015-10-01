# computing the probability of all the possible graphs or one specific graph 
prob = function( x, number.g = 4, G = NULL )
{
	if( !is.null( x $ phat ) ) stop( "Function needs output of 'bdgraph' with option save.all = TRUE" ) 
	
	sampleGraphs  = x $ sampleGraphs
	graphWeights  = x $ graphWeights
	sort_gWeights = sort( graphWeights, decreasing = TRUE )

	if( is.null(G) )
	{
		p      <- nrow( x $ lastGraph )
		list_G <- list()
		vec_G  <- c( rep( 0, p * ( p - 1 ) / 2 ) )  

		for ( i in 1 : number.g )
		{
			vec_G <- 0 * vec_G
			indG_i = sampleGraphs[ which( graphWeights == sort_gWeights[i] ) ]
			vec_G[ which( unlist( strsplit( as.character( indG_i ), "" ) ) == 1 ) ] <- 1
			list_G[[i]] <- matrix( 0, p, p )
			list_G[[i]][ upper.tri( list_G[[i]] ) ] <- vec_G
			list_G[[i]] <- Matrix( list_G[[i]], sparse = TRUE )
		}

		return( list( selected_G = list_G, prob_G = sort_gWeights[1 : number.g] / sum( graphWeights ) ) )
		
	} 
	else 
	{
		if ( class(G) == "sim" ) G <- G $ G

		G      = as.matrix(G)
		indG   = paste( G[upper.tri(G)], collapse = '' )
		wh     = which( sampleGraphs == indG )
		prob_G = ifelse( length(wh) == 0, 0, graphWeights[wh] / sum( graphWeights ) )

		return( prob_G )
	}
}
    
