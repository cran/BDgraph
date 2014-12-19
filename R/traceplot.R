# plot of graph size to check the convergency of BDMCMC algorithm
traceplot = function(output, acf = FALSE, pacf = FALSE, main = NULL, ...)
{
    allGraphs   <- output $ allGraphs
	graphWeights <- output $ graphWeights
	bestg   <- output $ sampleGraphs[which(max(graphWeights) == graphWeights)]	
	lin     <- length(which(unlist(strsplit(as.character(bestg), "")) == 1))
    y       <- sapply(allGraphs, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
	
	if (is.null(main)) main = "Trace of graph size"
	
	if (acf == FALSE & pacf == FALSE)
	{
		plot(x = 1 : length(allGraphs), y, type = "l", main = main,
			ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	   
	}
	
	if (acf == TRUE & pacf == TRUE)
	{
		op <- par(mfrow = c(2, 2), pty = "s") 
		plot(x = 1 : length(allGraphs), y, type = "l", main = main,
			   ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	  
		acf(y,  main = "ACF for graph size")
		pacf(y, main = "PACF for graph size")
		par(op)
	}
	
	if (acf == TRUE & pacf == FALSE)
	{
		op <- par(mfrow = c(1, 2), pty = "s") 
		plot(x = 1 : length(allGraphs), y, type = "l", main = main,
			   ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	  
		acf(y, main = "ACF for graph size")
		par(op)
	}
	
	if (acf == FALSE & pacf == TRUE)
	{
		op <- par(mfrow = c(1, 2), pty = "s") 
		plot(x = 1 : length(allGraphs), y, type = "l", main = main,
			   ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	  
		pacf(y, main = "PAIC for graph size")
		par(op)
	}		
}  
