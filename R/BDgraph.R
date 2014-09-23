## Main function: BDMCMC algorithm for selecting the best graphs 
# bdgraph with missing data option in copula
################################################################################
# in option "copulaNA", I should check ", PACKAGE = BDgraph)"
################################################################################
bdgraph = function( data, n = NULL, method = "exact", iter = 5000, 
					burnin = iter / 2, b = 3, D = NULL, Gstart = "empty" )
{
	startTime <- Sys.time()
	burnin = floor( burnin )
	
	if ( class(data) == "simulate" ) data <- data $ data

	if ( is.matrix(data) == FALSE & is.data.frame(data) == FALSE ) stop( "Data should be a matrix or dataframe" )
	if ( is.data.frame(data) ) data <- data.matrix(data)
	if ( any( is.na(data) ) )  method = "copulaNA"
	if ( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )

	dimd <- dim(data)
	p    <- dimd[2]
	n    <- dimd[1]

	if ( method == "copula" | method == "copulaNA" | method == "copula1" | method == "copulaNA1" )
	{
		Z              <- qnorm( apply( data, 2, rank, ties.method = "random" ) / (n + 1) )
		Zfill          <- matrix( rnorm( n * p ), n, p )   # for missing values
		Z[is.na(data)] <- Zfill[is.na(data)]               # for missing values
		Z              <- t( ( t(Z) - apply( Z, 2, mean ) ) / apply( Z, 2, sd ) )
		S              <- t(Z) %*% Z
		# ?? I should check it 
		R <- 0 * data
		for ( j in 1:p ) 
			R[,j] <- match( data[ , j], sort( unique( data[ , j] ) ) ) 
		R[is.na(R)] = 0 # dealing with missing values
	} else {
		if ( isSymmetric(data) )
		{
			if ( is.null(n) ) stop( "Please specify the number of observations 'n'" )
			cat( "The input is identified as the covriance matrix. \n" )
			S <- data
		} else {
			S <- t(data) %*% data
		}
	}

	if ( is.null(D) ) 
	{
		method = "exact1"
		if ( method == "copula" )   method == "copula1"
		if ( method == "copulaNA" ) method == "copulaNA1"
	}

	if ( is.null(D) )
	{ 
		D  <- diag(p)
		Ti <- D
		H  <- D
	} else {
		Ti <- chol( solve(D) )
		H  <- Ti / t( matrix( rep( diag(Ti) ,p ), p, p ) )
	}

	bstar <- b + n
	Ds    <- D + S
	invDs <- solve(Ds)
	Ts    <- chol(invDs)
	Hs    <- Ts / t( matrix( rep( diag(Ts) ,p ), p, p ) )	

	if ( class(Gstart) == "bdgraph" ) 
	{
		G <- Gstart $ lastGraph
		K <- Gstart $ lastK
	} 

	if ( class(Gstart) == "simulate" ) 
	{
		G <- as.matrix( Gstart $ G )
		K <- as.matrix( Gstart $ K )
	} 
	
	if ( class(Gstart) == "character" && Gstart == "empty"  )
	{
		G = 0 * S
		
		K <- matrix( 0, p, p )
		# rgwish ( double G[], double T[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ts), K = as.double(K), as.integer(bstar), 
					  as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 
	}
	
	if ( class(Gstart) == "character" && Gstart == "full" )
	{
		G       = matrix(1, p, p)
		diag(G) = 0
		# K       = rwishCpp( Ti = Ts, p = p, b = bstar )
		K = matrix( 0, p, p)
		# rwish ( double T[], double K[], int *p, int *b )
		result = .C( "rwish", as.double(Ts), K = as.double(K), as.integer(bstar), as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 
	}	

	if ( is.matrix(Gstart) )
	{
		G = as.matrix( Gstart )
		diag(G) = 0
		
		K <- matrix( 0, p, p )
		# rgwish ( double G[], double T[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ts), K = as.double(K), as.integer(bstar), 
					  as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 	
	}
		
	allGraphs    <- c( rep ( "a", iter ) ) # vector of numbers like "10100"
	allWeights   <- c( rep ( 0, iter ) )   # waiting time for every state		

	sampleGraphs <- c( rep ( "a", iter ) ) # vector of numbers like "10100" 
	graphWeights <- c( rep ( 0, iter ) ) # waiting time for every state
	
	sizeSampleG = 0

	rates     = 0 * K
	Ksum      = rates
	lastGraph = rates
	lastK     = rates

	mes <- paste( c(" ", iter," iteration is started.                    " ), collapse = "" )
	cat( mes, "\r" )
	#flush.console()

	if ( method == "exact" )
	{	
# bdmcmcExact( int *iter, int *burnin, double G[], double T[], double Ts[], double K[], int *p, 
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 double lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcExact", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(T), as.double(Ts), as.double(K), as.integer(p), 
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds)
				    , PACKAGE = "BDgraph" )
		################################################################################
	}

	if ( method == "exact1" )
	{
# bdmcmcApprox( int *iter, int *burnin, double G[], double Ts[], double K[], int *p, 
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 double lastGraph[], double lastK[],
#			 int *b, int *bstar, double Ds[] )
		result = .C( "bdmcmcExact1", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(Ts), as.double(K), as.integer(p), 
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(Ds)
				    , PACKAGE = "BDgraph" )
		################################################################################	
	}

		
	if ( method == "approx" )
	{
# bdmcmcApprox( int *iter, int *burnin, double G[], double T[], double Ts[], double K[], int *p, 
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 double lastGraph[], double lastK[],
#			 int *b, int *bstar, double Ti[], double Ds[] )
		result = .C( "bdmcmcApprox", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(T), as.double(Ts), as.double(K), as.integer(p), 
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(Ti), as.double(Ds)
				    , PACKAGE = "BDgraph" )
		################################################################################	
	}

	if ( method == "copula" )
	{
# bdmcmcCopula( int *iter, int *burnin, double G[], double Ti[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 double lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopula", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(Ti), as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds)
				    , PACKAGE = "BDgraph" )
				################################################################################
	}

	if ( method == "copulaNA1" )
	{
# bdmcmcCopula( int *iter, int *burnin, double G[], double Ti[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 double lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopulaNA1", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds) 
				    , PACKAGE = "BDgraph" )
		################################################################################
	}
	
	if ( method == "copula1" )
	{
#void bdmcmcCopula1( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopula1", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds)
				    , PACKAGE = "BDgraph" )
				################################################################################
	}

	if ( method == "copulaNA1" )
	{
# bdmcmcCopula( int *iter, int *burnin, double G[], double Ti[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 double lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopulaNA", as.integer(iter), as.integer(burnin), as.integer(G), 
		            as.double(Ti), as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds) 
				    , PACKAGE = "BDgraph" )
		################################################################################
	}

	Ksum         = matrix( result $ Ksum, p, p )
	allGraphs    = result $ allGraphs
	allWeights   = result $ allWeights
	sizeSampleG  = result $ sizeSampleG
	sampleGraphs = result $ sampleGraphs[1:sizeSampleG]
	graphWeights = result $ graphWeights[1:sizeSampleG]
	lastGraph    = matrix( result $ lastGraph, p, p )
	lastK        = matrix( result $ lastK, p, p )

#	mes <- paste( c(" ", iter," iteration done.                 " ), collapse = "" )
#	cat( mes, "\r" )
#	cat( "\n" )
#	flush.console()
	print( Sys.time() - startTime )  

	colnames( lastGraph ) = colnames( data )
	output <- list( sampleGraphs = sampleGraphs, graphWeights = graphWeights, Khat = Ksum / (iter - burnin), 
					  allGraphs = allGraphs, allWeights = allWeights, lastGraph = lastGraph, lastK = lastK )

	class( output ) <- "bdgraph"
	return( output )   
}
   
# computing probability of all links of the graph
phat = function(output, round = 3)
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

	dimlab <- dimnames( output $ lastGraph ) # lastG
	if ( is.null( dimlab ) )
	{
		dimlab <- as.character(1 : p)
		phat   <- matrix(0, p, p, dimnames = list(dimlab, dimlab))
	} else {
		phat   <- matrix(0, p, p, dimnames = dimlab)
	}

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

	dimlab   <- dimnames( output $ lastGraph )
	if (is.null(dimlab))
	{ 
		dimlab <- as.character(1 : p)
		graphi <- matrix(0, p, p, dimnames = list(dimlab, dimlab))
	} else {
		graphi <- matrix(0, p, p, dimnames = dimlab)
	}	

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
# To compare the result
roc = function ( G, est ) 
{
	G   <- as.matrix(G)     # G is the adjacency matrix of true graph 
	est <- as.matrix(est)   # est is the adjacency matrix of estimated graph 
	G[ lower.tri( G, diag = TRUE ) ]     <- 0
	est[ lower.tri( est, diag = TRUE ) ] <- 0
	p      <- nrow(G)
	
	tp.all <- ( G != 0 ) * ( est != 0 ) 
	fp.all <- ( G == 0 ) * ( est != 0 ) 
	fn.all <- ( G != 0 ) * ( est == 0 ) 
	tn.all <- ( G == 0 ) * ( est == 0 )
	
	tp     <- sum( tp.all )
	fp     <- sum( fp.all )
	fn     <- sum( fn.all )
	tn     <- sum( tn.all[ upper.tri( tn.all == 1 ) ] )
	
	# positive predictive value  
	Precision <- tp / ( tp + fp ) 
	
	if ( is.na(Precision) ) Precision <- 0
	# Precision is the probability that a randomly selected link is relevant
	Recall <- tp / ( tp + fn ) # also called TPR
	
	if ( is.na(Recall) ) Recall <- 0
	# Recall is the probability that a randomly selected relevant link 
	# is retrieved in a search 
	FPR <- fp / ( fp + tn ) # False positive rate
	
	if ( is.na(FPR) ) FPR <- 0
	Accuracy <- ( tp + tn ) / (tp + tn + fp + fn)
	
	if ( is.na(Accuracy) ) Accuracy <- 0
	# Specificity <- tn / (tn + fp) # or 1 - false positive rate 
	# Sensitivity <- tp / (tp + fn) # or true positive rate
	# F1score <- 2 * (Precision * Recall) / (Precision + Recall)
	F1score <- ( 2 * tp ) / ( 2 * tp + fp + fn )
	if ( is.na(F1score) ) F1score <- 0
	# harmonic mean of precision and recall, called F-measure or balanced F-score:
	roc.matrix <- matrix( c(tp, tn, fp, fn, Recall, FPR, Accuracy, F1score, Precision), 9, 1 )
	
	return( round( roc.matrix, 3 ) )
}
# function for ROC plot
outRoc = function( G, prob, cut )
{
	G[ lower.tri( G, diag = TRUE ) ]     <- 0
	prob[ lower.tri( prob, diag = TRUE ) ] <- 0
	p = nrow(prob)
	pp = p * ( p - 1 ) / 2
	sumEdges = sum(G)
	sumNoEdges = pp - sumEdges
	
	tp = c( rep( 0, cut + 1 ) )
	tp[1] = 1
	fp = c( rep( 0, cut + 1 ) )
	fp[1] = 1

	cutPoint = (0:cut) / cut
	
	for ( i in 2:cut )
	{
		# checking for cut pints
		estG = matrix( 0, p, p )
		estG[prob > cutPoint[i]] = 1

		tp.all <- ( G != 0 ) * ( estG != 0 ) 
		fp.all <- ( G == 0 ) * ( estG != 0 ) 	
		tp[i]  <- sum( tp.all ) / sumEdges
		fp[i]  <- sum( fp.all ) / sumNoEdges
	}
	
	return( list( tp = tp, fp = fp ) )
}
# To plot ROC curve
plotroc = function( G, prob, prob2 = NULL, cut = 20, smooth = FALSE )
{
    if ( class(G)     == "simulate" ) G     <- as.matrix( G $ G )
    if ( class(prob)  == "bdgraph" )  prob  <- as.matrix( phat( prob, round = 10 ) ) 
    
    output = outRoc( G = G, prob = prob, cut = cut )
    x      = output $ fp
    y      = output $ tp
 	
	if ( smooth == TRUE )
	{
		fit = smooth.spline( x = x, y = y )
		x   = c( 0, fit $ x )
		y   = c( 0, fit $ y )
	}
	
	par( mar = c( 3.8, 4.2, 1.8, 1 ) )
    plot( x = x, y = y, type = "l", col = "black", lty = 1, cex.lab = 1.6, cex.main = 2.5, cex.axis = 1.7,
         main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate", ylim = c(0,1) )
  
    if( !is.null(prob2) )
    {
        if ( class(prob2)  == "bdgraph" ) prob2 <- as.matrix( phat( prob2, round = 10 ) ) 
        output2 = outRoc( G = G, prob = prob2, cut = cut )
		x2      = output2 $ fp
		y2      = output2 $ tp

		if ( smooth == TRUE )
		{
			fit2 = smooth.spline( x = x2, y = y2 )
			x2   = c( 0, fit2 $ x )
			y2   = c( 0, fit2 $ y )
		}
		
        points( x = x2, y = y2, type = "l", col = "blue", lty = 2, lw = 2,
                main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate", ylim = c(0,1) )
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
# Data generator according to the graph structure
bdgraph.sim = function( n = 2, p = 10, graph = "random", size = NULL, prob = 0.2, 
                        class = NULL, type = "Gaussian", cut = 4, b = 3, D = diag(p), 
                        K = NULL, sigma = NULL, mean = 0, vis = FALSE )
{
    if ( is.matrix(K) )  graph <- "fixed"
    
    if ( type == "normal" )     type = "Gaussian"
    if ( type == "non-normal" ) type = "non-Gaussian"
    
    if ( is.matrix(graph) )
	{
		G     <- graph
	    graph <- "fixed"
    } 
	
	# build the graph structure
	if ( graph == "random" )
	{
		G <- matrix( 0, p, p )

		if ( is.null(size) )
		{
			if ( prob < 0 | prob > 1 ) stop("'prob' should be between zero and one")
			
			G[upper.tri(G)] <- rbinom( p * (p - 1) / 2, 1, prob )
			
		} else {
			if ( size < 0 | size > p * (p - 1) / 2 )  stop("Graph size should be between zero and p*(p-1)/2")
			
			smp <- sample( 1 : (p * (p - 1) / 2), size, replace = FALSE )
			G[upper.tri(G)][smp] <- 1
		}
		
		G <- G + t(G)
	}
	
	if ( graph == "cluster" )
	{
		# partition variables
		if ( is.null(class) )
		{ 
			if ( !is.null(size) )
			{
				class <- length(size)
			} else {
				class <- max( 2, ceiling(p / 20) )
			}
		}

		g.large <- p %% class
		g.small <- class - g.large
		n.small <- floor( p / class )
		n.large <- n.small + 1
		vp      <- c( rep( n.small, g.small ), rep( n.large, g.large ) )
		 
		G       <- matrix(0, p, p)
		 
		if ( is.null(size) )
		{
			if ( is.null(prob) ) prob <- 0.2
			if ( class > 1 ) prob <- class * prob
			if (prob < 0 | prob > 1) stop( "'prob' should be between zero and one" )

			for ( i in 1 : class )
			{
				tmp <- if ( i == 1 ) (1 : vp[1]) else ( ( sum(vp[1 : (i-1)]) + 1 ) : sum(vp[1:i]) )
				gg                <- matrix( 0, vp[i], vp[i] )
				gg[upper.tri(gg)] <- rbinom( vp[i] * (vp[i] - 1) / 2, 1, prob )
				G[tmp, tmp]       <- gg
			}
		} else {

			if ( class != length(size) )  stop( "Number of graph sizes is not match with number of clusters" )
			if ( sum(size) < 0 | sum(size) > p * (p - 1) / 2 )   stop( "Total graph sizes should be between zero and p*(p-1)/2" )

			for (i in 1 : class)
			{
				tmp <- if ( i == 1 ) ( 1 : vp[1] ) else ( (sum(vp[1 : (i-1)]) + 1) : sum(vp[1:i]) )
				gg  <- matrix(0, vp[i], vp[i])
				smp <- sample( 1 : (vp[i] * (vp[i] - 1) / 2), size[i], replace = FALSE )
				gg[upper.tri(gg)][smp] <- 1
				G[tmp, tmp]            <- gg
			}
		}
		
		G <- G + t(G)	   
	}

	if( graph == "hub" )
	{
		# partition variables
		if ( is.null(class) ) class <- ceiling(p / 20) 

		# partition variables into groups
		g.large <- p %% class
		g.small <- class - g.large
		n.small <- floor( p / class )
		n.large <- n.small + 1
		g.list  <- c( rep(n.small, g.small), rep(n.large, g.large) )
		g.ind   <- rep( c(1:class), g.list )
		
		G <- matrix(0, p, p)
		
		for ( i in 1:class )
		{
			tmp           <- which( g.ind == i )
			G[tmp[1],tmp] <- 1
			G[tmp,tmp[1]] <- 1
		}
	}
	
	if ( graph == "circle" )
	{
	    G       <- toeplitz( c( 0, 1, rep( 0, p - 2 ) ) )
        G[1, p] <- 1
		G[p, 1] <- 1
	}

	if ( graph == "scale-free" )
	{
#		data_sim <- huge.generator( n = 2, d = 5, graph = "scale-free" )
		G = matrix(0, p, p)
		# scaleFree( int *G, int *p )
		resultGraph = .C( "scaleFree", G = as.integer(G), as.integer(p), PACKAGE = "BDgraph" )
		G           = matrix( resultGraph $ G, p, p ) 
	}
	
    if ( !is.null(sigma) ) K <- solve(sigma)   

    if ( is.matrix(K) )
    { 
		G     <- 1 * ( abs(K) > 0.02 )
		if( is.null(sigma) ) sigma <- solve(K)	
    } else {
		Ti      <- chol( solve(D) )
		diag(G) <- 0
		
		K <- matrix( 0, p, p )
		# rgwish ( double G[], double T[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ti), K = as.double(K), as.integer(b), 
					  as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 		
			
		sigma <- solve(K)
	}
	
	diag(G) <- 0
	p       <- nrow(G)
	
	# generate multivariate normal data
	if ( typeof(mean) == "double" ) mean <- rep(mean, p)
	R <- chol(sigma)
	z <- matrix( rnorm(p * n), p, n )
	d <- t(R) %*% z + mean
	d <- t(d)

	if ( type == "mixed" )
	{
		# generating mixed data which are 'count', 'ordinal', 'non-Gaussian', 
		# 'binary', and 'Gaussian', respectively.
		ps = floor( p / 5 )
		
		# generating count data
		col_number     <- c( 1:ps )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qpois( p = prob, lambda = 10 )

		# generating ordinal data
		col_number     <- c( (ps + 1):(2 * ps) )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qpois( p = prob, lambda = 2 )

		# generating non-Guassian data
		col_number     <- c( (2 * ps + 1):(3 * ps) )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qexp( p = prob, rate = 10 )

		# for binary data
		col_number     <- c( (3 * ps + 1):(4 * ps) )
		prob           <- pnorm( d[, col_number] )
		d[,col_number] <- qbinom(p = prob, size = 1, prob = 0.5 )
	}

	if ( type == "non-Gaussian" )
	{
		# generating multivariate continuous non-Gaussian data  
		prob <- pnorm( d )
		d    <- qexp( p = prob, rate = 10 )
	}

	if ( type == "discrete" )
	{
		runif_m   <- matrix( runif(cut * p), nrow = p, ncol = cut )   
		marginals <- apply( runif_m, 1, function(x) { qnorm( cumsum( x / sum(x) )[-length(x)] ) } )
			 
		for ( j in 1:p )
		{
			breaks <- c( min(d[,j]) - 1, marginals[,j], max(d[,j]) + 1 )  
			d[,j]  <- as.integer( cut( d[,j], breaks = breaks, right = FALSE ) )
		}	
	}
	
	# graph visualization
	if ( vis == TRUE )
	{
		graphG <- graph.adjacency( G, mode = "undirected", diag = FALSE )
		
		if ( p < 20 )
		{
			plot.igraph( graphG, layout = layout.circle, main = "True graph structure",
						edge.color = 'black', vertex.color = "white", vertex.size = 10 )
		} else {
			plot.igraph( graphG, layout = layout.circle, main = "True graph structure", 
						edge.color = 'black', vertex.color = "white", vertex.size = 2 )
		}
	}
	
	dimlab     <- as.character(1 : p)
	simulation <- list( G = Matrix(G, sparse = TRUE, dimnames = list(dimlab, dimlab)), 
	                    data = d, sigma = sigma, K = K, graph = graph, type = type )
	
	class(simulation) <- "simulate"
	return(simulation)
}
# Print function for simulation data
print.simulate = function (x, ...)
{
	p <- ncol(x $ sigma)

	cat( paste( "  Data generated by bdgraph.sim"                           ), fill = TRUE )
	cat( paste( "  Data type       =", x $ type                             ), fill = TRUE )
	cat( paste( "  Sample size     =", nrow(x $ data)                       ), fill = TRUE )
	cat( paste( "  Graph type      =", x $ graph                            ), fill = TRUE )
	cat( paste( "  Number of nodes =", p                                    ), fill = TRUE )
	cat( paste( "  Graph size      =", sum(x $ G) / 2                       ), fill = TRUE )
	cat( paste( "  Sparsity        =", round(sum(x $ G) / (p * (p - 1)), 4) ), fill = TRUE )
}
# plot for class "simulate" from bdgraph.sim function
plot.simulate = function(x, main = NULL, layout = layout.circle, ...)
{
    if (is.null(main)) main <- "True graph structure"
  	g <- graph.adjacency( as.matrix(x $ G), mode = "undirected", diag = FALSE )
	
    plot.igraph(g, main = main, layout = layout, ...)
}		
# plot for class bdgraph
plot.bdgraph = function(x, g = 1, layout = layout.circle, ...)
{
	list.G  <- x $ sampleGraphs
	graphWeights <- x $ graphWeights
	p       <- nrow( x $ lastGraph )
	prob.G  <- graphWeights / sum(graphWeights)
	graphi  <- list()
	gv      <- c(rep(0, p * (p - 1) / 2))

	if (g == 2) op <- par(mfrow = c(1, 2), pty = "s")
	if (g > 2 & g < 7)  op <- par(mfrow = c(2, g %% 2 + trunc(g / 2)), pty = "s")

	for (i in 1 : g)
	{
		if (g > 6) dev.new()  
		gi <- list.G[which(prob.G == sort(prob.G, decreasing = T)[i])]
		gv <- 0 * gv
		gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
		graphi[[i]] <- matrix(0, p, p)
		graphi[[i]][upper.tri(graphi[[i]])] <- gv
		dimnames(graphi[[i]]) <- dimnames( x $ lastGraph )
		G    <- graph.adjacency(graphi[[i]], mode = "undirected", diag = FALSE)
		main <- ifelse (i == 1, "Graph with highest probability", paste(c(i, "th graph"), collapse = ""))
		plot.igraph(G, layout = layout, main = main, sub = paste(c("Posterior probability = ", 
		            round(sort(prob.G, decreasing = TRUE)[i], 4)), collapse = ""), ...)	   
	}
	
	if (g > 1 & g < 7) par(op)
}
# summary of bdgraph output
summary.bdgraph = function( object, vis = TRUE, layout = layout.circle, ... )
{
	sampleGraphs <- object $ sampleGraphs
	graphWeights <- object $ graphWeights
	p            <- nrow( object $ lastGraph )
	gv           <- c( rep( 0, p * (p - 1 ) / 2) )

	dimlab   <- dimnames( object $ lastGraph )
	if ( is.null(dimlab) )
	{ 
		dimlab <- as.character(1 : p)
		graphi <- matrix(0, p, p, dimnames = list(dimlab, dimlab))
	} else {
		graphi <- matrix(0, p, p, dimnames = dimlab)
	}	

	prob.G <- graphWeights / sum(graphWeights)
	gi     <- sampleGraphs[which(prob.G == max(prob.G))]
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
	graphi[upper.tri(graphi)] <- gv 

	if ( vis )
	{
		# plot best graph
		G  <- graph.adjacency( graphi, mode = "undirected", diag = FALSE )
		 
		op <- par(mfrow = c(2, 2), pty = "s")

		plot.igraph(G, layout = layout, main = "Best graph",
		  sub = paste(c("Posterior probability = ", round(max(prob.G), 4)), collapse = ""), ...)
		
		# plot posterior distribution of graph
		plot(x = 1 : length(graphWeights), y = graphWeights / sum(graphWeights), type = "h", main = "Posterior probability",
			 ylab = "Pr(graph|data)", xlab = "graph")
		abline(h = max(graphWeights) / sum(graphWeights), col = "red")
		text(which(max(graphWeights) == graphWeights), max(graphWeights) / sum(graphWeights), "P(best graph|data)", col = "gray60", adj = c(0, + 1))
		
		# plot posterior distribution of graph size
		suma     <- sapply( sampleGraphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )
		xx       <- unique(suma)
		weightsg <- vector()

		for (i in 1 : length(xx))
		{
			weightsg[i] <- sum( graphWeights[which( suma == xx[i] )] )
		}

		plot( x = xx, y = weightsg / sum(graphWeights), type = "h", main = "Posterior probability",
			 ylab = "Pr(graph size|data)", xlab = "graph size" )

		# plot trace of graph size
		allGraphs <- object $ allGraphs
		yy        <- sapply( allGraphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )
		
		plot( x = 1 : length(allGraphs), yy, type = "l", main = "Trace for graph size",
			  ylab = "graph size", xlab = "iteration")
		
		abline(h = sum(graphi), col = "red")	  
		
		par(op)
	}
	
	# phat
	pvec <- 0 * gv
	for (i in 1 : length(sampleGraphs))
	{
		inp       <- which(unlist(strsplit(as.character(sampleGraphs[i]), "")) == 1)
		pvec[inp] <- pvec[inp] + graphWeights[i]
	}

	phat                  <- 0 * graphi
	phat[upper.tri(phat)] <- pvec / sum(graphWeights)
	# estimation for precision matrix 
	Khat        <- object $ Khat

	return.list <- list(best.graph = Matrix(graphi + t(graphi), sparse = TRUE), phat = Matrix(round(phat, 2), 
					  sparse = TRUE), Khat = round(Khat, 3))
					  
	return(return.list)
}
# print of the bdgraph output
print.bdgraph = function(x, round = 3, Khat = FALSE, phat = FALSE, ...)
{
	sampleGraphs <- x $ sampleGraphs
	graphWeights  <- x $ graphWeights
	p        <- nrow( x $ lastGraph )
	# best graph
	prob.G   <- graphWeights / sum(graphWeights)
	gv       <- c(rep(0, p * (p - 1) / 2))
	gi       <- sampleGraphs[which(prob.G == max(prob.G))]
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1

	dimlab   <- dimnames( x $ lastGraph )
	if (is.null(dimlab))
	{ 
		dimlab <- as.character(1 : p)
		graphi <- matrix(0, p, p, dimnames = list(dimlab, dimlab))
	} else {
		graphi <- matrix(0, p, p, dimnames = dimlab)
	}	

	graphi[upper.tri(graphi)] <- gv
	cat(paste(""), fill = TRUE)
	cat(paste("Adjacency matrix of best graph"), fill = TRUE)
	cat(paste(""), fill = TRUE)
	printSpMatrix(Matrix(graphi + t(graphi), sparse = TRUE), col.names = TRUE, note.dropping.colnames = FALSE)

	cat(paste(""), fill = TRUE)
	cat(paste("Size of best graph =", sum(graphi)), fill = TRUE)
	cat(paste("Posterior probability of best graph = ", round(max(graphWeights) / sum(graphWeights), round)), fill = TRUE)  
	cat(paste(""), fill = TRUE)

	# print for precision matrix
	if (Khat == TRUE)
	{
		cat(paste(""), fill = TRUE)
		cat(paste("Estimation of precision matrix"), fill = TRUE)
		cat(paste(""), fill = TRUE)
		print(round(x $ Khat, round))
	}

	# print for phat
	if (phat == TRUE)
	{
		pvec <- 0 * gv
		for (i in 1 : length(sampleGraphs))
		{
			inp       <- which(unlist(strsplit(as.character(sampleGraphs[i]), "")) == 1)
			pvec[inp] <- pvec[inp] + graphWeights[i]
		}

		phat                  <- 0 * graphi
		phat[upper.tri(phat)] <- pvec / sum(graphWeights)
		cat(paste(""), fill = TRUE)
		cat(paste("Posterior probability of links"), fill = TRUE)
		cat(paste(""), fill = TRUE)

		printSpMatrix(Matrix(round(phat, round), sparse = TRUE), col.names = TRUE, note.dropping.colnames = FALSE)  
	}
} 
# sampling from G-Wishart distribution
rgwish = function( n = 1, G = NULL, b = 3, D = NULL )
{
	if ( is.null(G) ) stop( "Adjacency matrix G should be determined" )
	G <- as.matrix(G)
	if ( sum( (G == 1) * (G == 0) ) != 0 ) stop( "Elements of matrix G should be zero or one" )
	
	G[lower.tri(G, diag(TRUE))] <- 0
	p <- nrow(G)  
	
	if ( is.null(D) ) D <- diag(p)	

	samples <- array( 0, c( p, p, n ) )

	for ( i in 1 : n )
	{
		Ti           <- chol( solve(D) ) 
		# samples[,,i] <- rgwish.exact(G = G + t(G), b = b, T = Ti, p = p, threshold = 1e-8)
		# samples[,,i] <- rgwishCpp( G = G + t(G), Ti = Ti, p = p, b = b )
		K = matrix( 0, p, p )
		G = G + t(G)
		# rgwish ( double G[], double T[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ti), K = as.double(K), as.integer(b), 
					  as.integer(p), PACKAGE = "BDgraph" )
		samples[,,i] = matrix ( result $ K, p, p ) 		
	}	

	return( samples )   
}
# sampling from Wishart distribution
rwish = function( n = 1, p = 2, b = 3, D = diag(p) )
{
	samples <- array( 0, c( p, p, n ) )

	for ( i in 1 : n )
	{
		Ti = chol( solve(D) ) 
		
		K  = matrix( 0, p, p )
		# rwish ( double T[], double K[], int *p, int *b )
		result = .C( "rwish", as.double(Ti), K = as.double(K), as.integer(b), as.integer(p) )
		
		samples[,,i] = matrix ( result $ K, p, p ) 		
	}	

	return( samples )   
}
# non-parametric transfer function for non-normal data
bdgraph.npn = function(data, npn = "shrinkage", npn.thresh = NULL)
{
    if (is.matrix(data) == FALSE & is.data.frame(data) == FALSE) stop("Data should be a matrix or dataframe")
	
    if (is.data.frame(data) == TRUE) data <- data.matrix(data)
	
    if (any(is.na(data))) stop("Data should contain no missing data") 
	
	n <- nrow(data)
  	# shrinkage transfer
	if(npn == "shrinkage")
	{
		data <- qnorm(apply(data, 2, rank) / (n + 1))
		data <- data / sd(data[ , 1])
	}
	
	# truncation transfer
	if(npn == "truncation")
	{
		if(is.null(npn.thresh)) npn.thresh <- 0.25 * (n ^ - 0.25) * (pi * log(n)) ^ - 0.5
		data <- qnorm( pmin(pmax(apply(data, 2, rank) / n, npn.thresh), 1 - npn.thresh) )
    	data <- data / sd(data[ , 1])
	}

	if(npn == "skeptic") data <- 2 * sin( pi / 6 * cor(data, method = "spearman") )
	
	return(data)
}
# Monte Carlo approximation of expectation in normalizing constant
log.Exp.MC = function(A, b, H, mc, p)
{
	nu  <- rowSums(A)  
	f_T <- c(rep(0, mc))

	for (k in 1 : mc)
	{
		psi         <- diag(sqrt(rchisq(p, b + nu)))
		psi[A == 1] <- rnorm(1)

		if (identical(H, diag(p)))
		{
			for (i in 2 : (p - 1))
			{
				for (j in (i + 1) : p)
				{
					if (A[i, j] == 0)
					{
						psi[i, j] <- - sum(psi[1 : (i - 1), i] * psi[1 : (i - 1), j]) / psi[i, i]
						f_T[k]    <- f_T[k] + psi[i, j] ^ 2
					}
				}
			}
			
		} else {
			for (i in 1 : (p - 1))
			{
				for (j in (i + 1) : p)
				{
					if (A[i, j] == 0)
					{
						psi[i, j] <- - sum(psi[i, i : (j - 1)] * H[i : (j - 1), j])

						if (i > 1)
						{
							for (r in 1 : (i - 1))
							{
								psi[i, j] <- psi[i, j] - ((sum(psi[r, r : i] * H[r : i, i])) *
									   (sum(psi[r, r : j] * H[r : j, j]))) / (psi[i, i])
							}
						}

						f_T[k] <- f_T[k] + psi[i, j] ^ 2
					}
				}
			}
		}
	}
	
	# check for infinity values
	f_T[!is.finite(f_T)] <- gamma(171)
	
	log.Exp <- log( mean( exp(- f_T / 2) ) )

	return( log.Exp )
}
# To compute Normalizing constant of G-Wishart distribution according to ATAY-KAYIS AND MASSAM (2005)
I.g = function( G, b = 3, D = diag( ncol(G) ), mc = 100 )
{
	if (b <= 2) stop("In G-Wishart distribution parameter 'b' has to be more than 2")

	G[lower.tri(G, diag = TRUE)] <- 0

	sumrowG <- rowSums(G)
	sumcolG <- colSums(G)
	p       <- nrow(G)

	Ti          <- chol(solve(D))
	H           <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
	log.Exp.f_T <- log.Exp.MC(G, b, H, mc, p)

	sumG    <- sum(G)
	c_dT    <- (sumG / 2) * log(pi) + (p * b / 2 + sumG) * log(2) +
	           sum(lgamma((b + sumrowG) / 2)) + sum((b + sumrowG + sumcolG) * log(diag(Ti)))
	  
	Ig      <- exp( c_dT + log.Exp.f_T )

	if ( is.finite(Ig) )
	{
		return( Ig )
		
	} else {
		logIg <- c_dT + log.Exp.f_T
		cat(paste(""), fill = TRUE)
		cat(paste("Normalizing constant is infinte"), fill = TRUE)
		cat(paste("Log of normalizing constant =", logIg), fill = TRUE)
		cat(paste(""), fill = TRUE)

		return( logIg ) 
	}
}






