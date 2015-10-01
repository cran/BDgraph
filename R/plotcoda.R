# To check the convergency of the BDMCMC algorithm
plotcoda = function( x, thin = NULL, main = NULL, links = TRUE, ... )
{
	if( !is.null( x $ phat ) ) stop( "Function needs output of 'bdgraph' with option save.all = TRUE" ) 
	
	if( is.null( thin ) ) thin = ceiling( length( x $ allGraphs ) / 1000 )

	sampleGraphs    = x $ sampleGraphs
	p               = nrow( x $ lastGraph )
	qp              = p * ( p - 1 ) / 2 
	allWeights      = x $ allWeights
	allGraphs       = x $ allGraphs

	allG_new        = allGraphs[ c( thin * ( 1 : floor( length( allGraphs ) / thin ) ) ) ]
	allWeights_new  = allWeights[ c( thin * ( 1 : floor( length( allWeights ) / thin ) ) ) ]
	length_allG_new = length( allG_new )
	result          = matrix( 0, qp, length_allG_new )
	vec_result      = 0 * result[ , 1]

	for ( g in 1 : length_allG_new )
	{
		mes = paste( c( "Calculation ... in progress : ", floor( 100 * g / length_allG_new ), "%" ), collapse = "" )
		cat(mes, "\r")
		flush.console()	

		which_edge             = which( unlist( strsplit( as.character( sampleGraphs[ allG_new[g] ] ), "" ) ) == 1 )
		vec_result[which_edge] = vec_result[which_edge] + allWeights_new[g]
		result[ ,g]            = vec_result / sum( allWeights_new[ c( 1 : g ) ] )    	 
	}

	if ( links )
		if ( p > 15 )
		{
			randomLinks = sample( x = 1:qp, size = ( qp - 100 ), replace = FALSE )
			result[ randomLinks, ] = 0
		}
	
	mes = paste( c( "Calculation ... done.                        " ), collapse = "" )
	cat(mes, "\r")
	cat("\n")
	flush.console()

	matplot( x = thin * ( 1 : length_allG_new ), y = t( result ), type = "l", lty = 1, col = 1,
		  xlab = "Iteration", ylab = "Posterior link probability", cex.lab = 1.3, cex.axis = 1.2 )
		  
	if ( is.null( main ) ) main = "Trace of the Posterior Probabilities of the Links"
	title( main = main, cex.main = 1.5 )
}
    
