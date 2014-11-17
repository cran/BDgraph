# computing probability of all links of the graph
phat = function( output, round = 3 )
{
	sampleGraphs <- output $ sampleGraphs
	graphWeights <- output $ graphWeights
	p            <- nrow( output $ lastGraph )
	pvec         <- c(rep(0, p * (p - 1) / 2))

	for (i in 1 : length(sampleGraphs))
	{
		inp       <- which( unlist( strsplit( as.character(sampleGraphs[i]), "" ) ) == 1 )
		pvec[inp] <- pvec[inp] + graphWeights[i]
	}

	dimlab <- colnames( output $ lastGraph ) # lastG
	if ( is.null( dimlab ) ) dimlab <- as.character(1 : p)
	
	phat   <- matrix( 0, p, p, dimnames = list(dimlab, dimlab) )
	phat[upper.tri(phat)] <- pvec / sum(graphWeights)

	return(Matrix(round(phat, round)))
}
# To check the convergency of the BDMCMC algorithm
plotcoda = function(output, thin = NULL, trace = TRUE, main = NULL, ...)
{
	if (is.null(thin)) thin = ceiling(length(output $ allGraphs) / 1000)

	op          <- par(mfrow = c(2, 2), pty = "s")
	p           <- nrow( output $ lastGraph ) 
	allWeights <- output $ allWeights
	allGraphs       <- output $ allGraphs

	graphWeights <- output $ graphWeights
	bestg   <- output $ sampleGraphs[which(max(graphWeights) == graphWeights)]	
	lin     <- length(which(unlist(strsplit(as.character(bestg), "")) == 1))
	y       <- sapply(allGraphs, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))

	plot(x = (1 : length(allGraphs)), y, type = "l", main = "Trace of graph size",
	   ylab = "graph size", xlab = "iteration", ...)
	abline(h = lin, col = "red")
	acf(y, main = "ACF for graph size")
	pacf(y, main = "PACF for graph size")

	allG.new        <- allGraphs[c(thin * (1 : floor(length(allGraphs) / thin)))]
	allWeights.new <- allWeights[c(thin * (1 : floor(length(allWeights) / thin)))]
	length.allG.new <- length(allG.new)
	ff              <- matrix(0, p * (p - 1) / 2, length.allG.new)
	ffv             <- 0 * ff[ , 1]

	for (g in 1 : length.allG.new)
	{
		if (trace == TRUE)
		{
			mes <- paste(c("Calculation ... in progress : ", floor(100 * g / length.allG.new), "%"), collapse = "")
			cat(mes, "\r")
			flush.console()	
		}

		inp      <- which(unlist(strsplit(as.character(allG.new[g]), "")) == 1)
		ffv[inp] <- ffv[inp] + allWeights.new[g]
		ff[ ,g]  <- ffv / sum(allWeights.new[c(1 : g)])    	 
	}

	if(trace == TRUE)
	{
		mes <- paste(c("Calculation ... done.                        "), collapse = "")
		cat(mes, "\r")
		cat("\n")
		flush.console()
	} 

	matplot(x = thin * (1 : length.allG.new), y = t(ff), type = "l", lty = 1, col = 1,
		  xlab = "iteration", ylab = "posterior link probability")
		  
	if (is.null(main)) main <- "Trace plot"
	title(main = main)
	#abline(v = thin * length.allG.new / 2, col = "blue")

	par(op)
}
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
# To select the best graph (graph with the highest posterior probability) 
select = function ( output, vis = FALSE )
{
	sampleGraphs   <- output $ sampleGraphs
	graphWeights    <- output $ graphWeights
	p          <- nrow( output $ lastGraph )
	prob.G     <- graphWeights / sum(graphWeights)
	max.prob.G <- which(prob.G == max(prob.G))
	
	if(length(max.prob.G) > 1) max.prob.G <- max.prob.G[1] 
	
	gi        <- sampleGraphs[max.prob.G]
	gv        <- c(rep(0, p * (p - 1) / 2))
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1

	dimlab   <- colnames( output $ lastGraph )
	if ( is.null(dimlab) ) dimlab <- as.character(1 : p)
	
	graphi <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	

	graphi[upper.tri(graphi)] <- gv
	if (vis == TRUE)
	{
		G <- graph.adjacency(graphi, mode = "undirected", diag = FALSE)
		
		plot.igraph(G, layout = layout.circle, main = "Graph with highest probability", 
		    sub = paste(c("Posterior probability = ", round(max(prob.G), 4)), collapse = ""))
	}

	return(Matrix(graphi + t(graphi), sparse = TRUE))
}
# computing the probability of all the possible graphs or one specific graph 
prob = function( output, g = 4, G = NULL )
{
	sampleGraphs <- output $ sampleGraphs
	graphWeights  <- output $ graphWeights

	if (is.null(G))
	{
		p      <- nrow( output $ lastGraph )
		graphi <- list()
		gv     <- c(rep(0, p * (p - 1) / 2))  

		for (i in 1 : g)
		{
			gi <- sampleGraphs[which(graphWeights == sort(graphWeights, decreasing = T)[i])]
			gv <- 0 * gv
			gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
			graphi[[i]] <- matrix(0, p, p)
			graphi[[i]][upper.tri(graphi[[i]])] <- gv
			graphi[[i]] <- Matrix(graphi[[i]] + t(graphi[[i]]), sparse = TRUE)
		}

		return(list(best.G = graphi, prob.G = sort(graphWeights, decreasing = T)[1 : g] / sum(graphWeights)))
		
	} else {
		if (class(G) == "simulate") G <- G $ G

		G     <- as.matrix(G)
		indA  <- paste(G[upper.tri(G)], collapse = '')
		wh    <- which(sampleGraphs == indA)
		probG <- ifelse(length(wh) == 0, 0, graphWeights[wh] / sum(graphWeights))

		return(probG)
	}
}
# To compare the result according to the true graph
compare = function ( G, est, est2 = NULL, colnames = NULL, vis = FALSE ) 
{
	if ( class(G)    == "simulate" ) G    <- G $ G 
	if ( class(G)    == "bdgraph" ){ es <- select(G) ; G <- est $ G ; est <- es }
	if ( class(est)  == "bdgraph" )  est  <- select( est ) 
	if ( class(est2) == "bdgraph" )  est2 <- select( est2 ) 
	if ( class(est)  == "select" )   est  <- est $ refit
	if ( class(est2) == "select" )   est2 <- est2 $ refit

	compare.true <- roc( G = G, est = G )
	compare.g    <- roc( G = G, est = est )

	if ( !is.null(est2) )
	{
		compare.g2  <- roc( G = G, est = est2 )
		compare.all <- cbind( compare.true, compare.g, compare.g2 )
		if ( is.null(colnames) ) colnames <- c("True graph", "estimate", "estimate2")
	} else {
		compare.all <- cbind( compare.true, compare.g )
		if ( is.null(colnames) ) colnames <- c( "True graph", "estimate" )
	}

   if ( vis == TRUE )
   {
		G   <- graph.adjacency( G,   mode = "undirected", diag = FALSE )
		est <- graph.adjacency( est, mode = "undirected", diag = FALSE )

		if ( is.null(est2) )
		{
			op <- par( mfrow = c(1, 2), pty = "s" )
			plot.igraph( as.matrix(G), layout = layout.circle, main = colnames[1] )
			plot.igraph( as.matrix(est), layout = layout.circle, main = colnames[2] )
		} else {
			op   <- par( mfrow = c(2, 2), pty = "s" )
			est2 <- graph.adjacency(as.matrix(est2), mode = "undirected", diag = FALSE)
			plot.igraph( G,    layout = layout.circle, main = colnames[1] )
			plot.igraph( est,  layout = layout.circle, main = colnames[2] )
			plot.igraph( est2, layout = layout.circle, main = colnames[3] )			
		}
		par(op)
   }
   
   colnames( compare.all ) <- colnames
   rownames( compare.all ) <- c("true positive", "true negative", "false positive", "false negative", 
                "true positive rate", "false positive rate", "accuracy", "balanced F-score", "positive predictive")
				
   return( compare.all )
}
