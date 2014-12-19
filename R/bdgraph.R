## Main function: BDMCMC algorithm for selecting the best graphs 
# bdgraph with missing data option in copula
################################################################################
# in option "copulaNA", I should check ", PACKAGE = BDgraph)"
################################################################################
bdgraph = function( data, n = NULL, method = "ggm", iter = 5000, 
					burnin = iter / 2, b = 3, D = NULL, Gstart = "empty" )
{
	startTime <- Sys.time()
	burnin = floor( burnin )
	
	if ( class(data) == "simulate" ) data <- data $ data

	if ( is.matrix(data) == FALSE & is.data.frame(data) == FALSE ) stop( "Data should be a matrix or dataframe" )
	if ( is.data.frame(data) ) data <- data.matrix(data)
	if ( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )
	if ( any( is.na(data) ) ) 
	{
		if ( method == "ggm" ) 
		{
			stop( "ggm method does not deal with missing value. You could choose method = cggm" )
		} else {
			method = "cggm-NA"
		}
	}
		
	dimd <- dim(data)
	p    <- dimd[2]
	n    <- dimd[1]

	if ( method == "cggm" | method == "cggm-NA" | method == "cggm-dmh" | method == "cggm-dmh-NA" )
	{
		Z              <- qnorm( apply( data, 2, rank, ties.method = "random" ) / (n + 1) )
		Zfill          <- matrix( rnorm( n * p ), n, p )   # for missing values
		Z[is.na(data)] <- Zfill[is.na(data)]               # for missing values
		Z              <- t( ( t(Z) - apply( Z, 2, mean ) ) / apply( Z, 2, sd ) )
		S              <- t(Z) %*% Z
		# ?? I should check
		R <- 0 * data
		for ( j in 1:p ) 
			R[,j] = match( data[ , j], sort( unique( data[ , j] ) ) ) 
		R[ is.na(R) ] = 0     # dealing with missing values
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
		D  <- diag(p)
		Ti <- D
	} else {
		Ti <- chol( solve(D) )
		if ( method == "ggm" )     method == "ggm-approx"
		if ( method == "cggm" )    method == "cggm-dmh"
		if ( method == "cggm-NA" ) method == "cggm-dmh-NA"
	}

	bstar <- b + n
	Ds    <- D + S
	invDs <- solve(Ds)
	Ts    <- chol(invDs)

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
		# rgwish ( integer G[], double Ti[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ts), K = as.double(K), as.integer(bstar), 
					  as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 
	}
	
	if ( class(Gstart) == "character" && Gstart == "full" )
	{
		G       = matrix(1, p, p)
		diag(G) = 0
	
		K = matrix( 0, p, p)
		# rwish ( double Ti[], double K[], int *p, int *b )
		result = .C( "rwish", as.double(Ts), K = as.double(K), as.integer(bstar), as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 
	}	

	if ( is.matrix(Gstart) )
	{
		G       = Gstart
		diag(G) = 0
		
		K <- matrix( 0, p, p )
		# rgwish ( integer G[], double Ti[], double K[], int *b, int *p )
		result = .C( "rgwish", as.integer(G), as.double(Ts), K = as.double(K), as.integer(bstar), 
					  as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 	
	}
		
	allGraphs    <- c( rep ( "a", iter ) ) # vector of numbers like "10100"
	allWeights   <- c( rep ( 0, iter ) )   # waiting time for every state		

	sampleGraphs <- c( rep ( "a", iter ) ) # vector of numbers like "10100" 
	graphWeights <- c( rep ( 0, iter ) )   # waiting time for every state
	
	sizeSampleG = 0

	rates     = 0 * K
	Ksum      = rates
	lastGraph = rates
	lastK     = rates

	mes <- paste( c(iter," iteration is started.                    " ), collapse = "" )
	cat( mes, "\r" )
	#flush.console()

	if ( method == "ggm" )
	{
# bdmcmcExact( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double Ds[] )   
		result = .C( "bdmcmcExact", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ts), as.double(K), as.integer(p), 
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(Ds)
				    , PACKAGE = "BDgraph" )
	}
	
	if ( method == "ggm-dmh" )
	{	
# bdmcmcDmh( int *iter, int *burnin, int G[], double Ti[], double Ts[], double K[], int *p, 
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcDmh", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ti), as.double(Ts), as.double(K), as.integer(p), 
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds)
				    , PACKAGE = "BDgraph" )
	}

	if ( method == "ggm-approx" )
	{
# bdmcmcApprox( int *iter, int *burnin, int G[], double Ti[], double Ts[], double K[], int *p, 
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double Ds[] )
		result = .C( "bdmcmcApprox", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ti), as.double(Ts), as.double(K), as.integer(p), 
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(Ds)
				    , PACKAGE = "BDgraph" )
	}

	if ( method == "cggm" )
	{
# bdmcmcCopula( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopula", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds)
				    , PACKAGE = "BDgraph" )
	}

	if ( method == "cggm-NA" )
	{
# bdmcmcCopulaNA( int *iter, int *burnin, int G[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopulaNA", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds) 
				    , PACKAGE = "BDgraph" )
	}

	if ( method == "cggm-dmh" )
	{
# bdmcmcCopulaDmh( int *iter, int *burnin, int G[], double Ti[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopulaDmh", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ti), as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds)
				    , PACKAGE = "BDgraph" )
	}

	if ( method == "cggm-dmh-NA" )
	{
# bdmcmcCopulaDmhNA( int *iter, int *burnin, int G[], double Ti[], double Ts[], double K[], int *p, 
#			 double Z[], int R[], int *n,
#			 string allGraphs[], double allWeights[], double Ksum[], 
#			 string sampleGraphs[], double graphWeights[], int *sizeSampleG,
#			 int lastGraph[], double lastK[],
#			 int *b, int *bstar, double D[], double Ds[] )
		result = .C( "bdmcmcCopulaDmhNA", as.integer(iter), as.integer(burnin), as.integer(G), as.double(Ti), as.double(Ts), as.double(K), as.integer(p),
		            as.double(Z), as.integer(R), as.integer(n),
					allGraphs = as.character(allGraphs), allWeights = as.double(allWeights), Ksum = as.double(Ksum), 
				    sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
				    lastGraph = as.integer(lastGraph), lastK = as.double(lastK),
				    as.integer(b), as.integer(bstar), as.double(D), as.double(Ds) 
				    , PACKAGE = "BDgraph" )
	}
	
	Ksum         = matrix( result $ Ksum, p, p )
	allGraphs    = result $ allGraphs
	allWeights   = result $ allWeights
	sizeSampleG  = result $ sizeSampleG
	sampleGraphs = result $ sampleGraphs[1:sizeSampleG]
	graphWeights = result $ graphWeights[1:sizeSampleG]
	lastGraph    = matrix( result $ lastGraph, p, p )
	lastK        = matrix( result $ lastK, p, p )

	print( Sys.time() - startTime )  

	colnames( lastGraph ) = colnames( data )
	output <- list( sampleGraphs = sampleGraphs, graphWeights = graphWeights, Khat = Ksum / (iter - burnin), 
					  allGraphs = allGraphs, allWeights = allWeights, lastGraph = lastGraph, lastK = lastK )

	class( output ) <- "bdgraph"
	return( output )   
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
		colnames(graphi[[i]]) = colnames( x $ lastGraph )
		rownames(graphi[[i]]) = colnames(graphi[[i]])
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

	dimlab   <- colnames( object $ lastGraph )
	if ( is.null(dimlab) ) dimlab <- as.character(1 : p)
	
	graphi <- matrix( 0, p, p, dimnames = list(dimlab, dimlab) )	

	prob.G <- graphWeights / sum(graphWeights)
	gi     <- sampleGraphs[which(prob.G == max(prob.G))]
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
	graphi[upper.tri(graphi)] <- gv 

	if ( vis )
	{
		# plot best graph
		G  <- graph.adjacency( graphi, mode = "undirected", diag = FALSE )
		 
		op = par( mfrow = c(2, 2), pty = "s", omi = c(0.3,0.3,0.3,0.3), mai = c(0.3,0.3,0.3,0.3) ) 

		plot.igraph(G, layout = layout, main = "Best graph",
		  sub = paste( c( "Posterior probability = ", max( prob.G ) ), collapse = "" ), ... )
		
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

	dimlab   <- colnames( x $ lastGraph )
	if ( is.null(dimlab) ) dimlab <- as.character(1 : p)
	
	graphi <- matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	

	graphi[upper.tri(graphi)] <- gv
	cat(paste(""), fill = TRUE)
	cat(paste("Adjacency matrix of best graph"), fill = TRUE)
	cat(paste(""), fill = TRUE)
	printSpMatrix(Matrix(graphi + t(graphi), sparse = TRUE), col.names = TRUE, note.dropping.colnames = FALSE)

	cat( paste( "" ), fill = TRUE )
	cat( paste( "Size of best graph =", sum(graphi) ), fill = TRUE )
	cat( paste( "Posterior probability of best graph = ", max(graphWeights) / sum(graphWeights) ), fill = TRUE )  
	cat( paste( ""), fill = TRUE )

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
   





