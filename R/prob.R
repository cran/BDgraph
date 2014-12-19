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
