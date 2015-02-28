# To check the convergency of the BDMCMC algorithm
plotcoda = function( output, thin = NULL, main = NULL, ... )
{
	if ( is.null(thin) ) thin = ceiling( length( output $ allGraphs ) / 1000 )

	p          <- nrow( output $ lastGraph ) 
	allWeights <- output $ allWeights
	allGraphs  <- output $ allGraphs

	allG.new        <- allGraphs[c(thin * (1 : floor(length(allGraphs) / thin)))]
	allWeights.new <- allWeights[c(thin * (1 : floor(length(allWeights) / thin)))]
	length.allG.new <- length(allG.new)
	ff              <- matrix(0, p * (p - 1) / 2, length.allG.new)
	ffv             <- 0 * ff[ , 1]

	for ( g in 1 : length.allG.new )
	{
		mes <- paste( c( "Calculation ... in progress : ", floor( 100 * g / length.allG.new ), "%" ), collapse = "" )
		cat(mes, "\r")
		flush.console()	

		inp      <- which( unlist(strsplit(as.character(allG.new[g]), "")) == 1 )
		ffv[inp] <- ffv[inp] + allWeights.new[g]
		ff[ ,g]  <- ffv / sum(allWeights.new[c(1 : g)])    	 
	}

	mes <- paste(c("Calculation ... done.                        "), collapse = "")
	cat(mes, "\r")
	cat("\n")
	flush.console()

	matplot( x = thin * (1 : length.allG.new), y = t(ff), type = "l", lty = 1, col = 1,
		  xlab = "Iteration", ylab = "Posterior link probability", cex.lab = 1.3, cex.axis = 1.2 )
		  
	if ( is.null(main) ) main <- "Trace of the Posterior Probabilities of the Links"
	title( main = main, cex.main = 1.5 )
}
    
