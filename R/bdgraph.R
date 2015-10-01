## Main function: BDMCMC algorithm for graphical models 
################################################################################
bdgraph = function( data, n = NULL, method = "ggm", algorithm = "bdmcmc", 
					iter = 5000, burnin = iter / 2, b = 3, Gstart = "empty",
					save.all = FALSE )
{
	threshold = 1e-8  # for sampling from gwishart distribution
	startTime = Sys.time()
	burnin    = floor( burnin )
	
	if( class(data) == "sim" ) data <- data $ data

	if( !is.matrix(data) & !is.data.frame(data) ) stop( "Data should be a matrix or dataframe" )
	if( is.data.frame(data) ) data <- data.matrix(data)
	if( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )

	if( any( is.na(data) ) ) 
	{
		if( method == "ggm" ) stop( "ggm method does not deal with missing value. You could choose method = gcgm" )	
		gcgm_NA = 1
	}
	else
	{
		gcgm_NA = 0
	}
		
	dimd <- dim(data)
	p    <- dimd[2]
	if( is.null(n) ) n <- dimd[1]

	if( method == "gcgm" )
	{
		R <- 0 * data
		for( j in 1:p ) R[,j] = match( data[ , j], sort( unique( data[ , j] ) ) ) 
		R[ is.na(R) ] = 0     # dealing with missing values	

		# copula for continuous non-Gaussian data
		if( gcgm_NA == 0 && min( apply( R, 2, max ) ) > ( n - 5 * n / 100 ) )
		{
			# copula transfer 
			data = qnorm( apply( data, 2, rank ) / ( n + 1 ) )
#~ 			data = data / sd( data[ , 1] )
		
			method = "ggm"
		}
		else
		{	# for non-Gaussian data
			Z              <- qnorm( apply( data, 2, rank, ties.method = "random" ) / (n + 1) )
			Zfill          <- matrix( rnorm( n * p ), n, p )   # for missing values
			Z[is.na(data)] <- Zfill[is.na(data)]               # for missing values
			Z              <- t( ( t(Z) - apply( Z, 2, mean ) ) / apply( Z, 2, sd ) )
			S              <- t(Z) %*% Z
		}
	} 
	
	if( method == "ggm" ) 
	{
		if( isSymmetric(data) )
		{
			if ( is.null(n) ) stop( "Please specify the number of observations 'n'" )
			cat( "Input is identified as the covriance matrix. \n" )
			S <- data
		}
		else
		{
 			S <- t(data) %*% data
		}
	}

	D     = diag(p)
	bstar <- b + n
	Ds    <- D + S
	Ts    <- chol( solve( Ds ) )

	if( class(Gstart) == "bdgraph" ) 
	{
		G <- Gstart $ lastGraph
		K <- Gstart $ lastK
	} 

	if( class(Gstart) == "sim" ) 
	{
		G <- as.matrix( Gstart $ G )
		K <- as.matrix( Gstart $ K )
	} 
	
	if( class(Gstart) == "character" && Gstart == "empty"  )
	{
		G = 0 * S
		K = G
		
		result = .C( "rgwish", as.integer(G), as.double(Ts), K = as.double(K), as.integer(bstar), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 
	}
	
	if( class(Gstart) == "character" && Gstart == "full" )
	{
		G       = matrix(1, p, p)
		diag(G) = 0
		K       = 0 * G

		result = .C( "rwish", as.double(Ts), K = as.double(K), as.integer(bstar), as.integer(p), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 
	}	

	if( is.matrix( Gstart ) )
	{
		G       = Gstart
		diag(G) = 0
		K       = 0 * G
		
		result = .C( "rgwish", as.integer(G), as.double(Ts), K = as.double(K), as.integer(bstar), as.integer(p), as.double(threshold), PACKAGE = "BDgraph" )
		K = matrix ( result $ K, p, p ) 	
	}
		
	if( save.all == TRUE )
	{
		qp1          = ( p * ( p - 1 ) / 2 ) + 1
		stringG      = paste( c( rep( 0, qp1 ) ), collapse = '' )
		sampleGraphs = c( rep ( stringG, iter - burnin ) )  # vector of numbers like "10100" 
		graphWeights = c( rep ( 0, iter - burnin ) )        # waiting time for every state
		allGraphs    = c( rep ( 0, iter - burnin ) )        # vector of numbers like "10100"
		allWeights   = c( rep ( 1, iter - burnin ) )        # waiting time for every state		
		sizeSampleG  = 0
	}
	else
	{
		phat = 0 * K
	}

	if( ( save.all == TRUE ) && ( p > 50 & iter > 20000 ) )
	{
		cat( "  WARNING: Memory needs to run this function is around " )
		print( ( iter - burnin ) * object.size( stringG ), units = "auto" ) 
	} 
	
	Khat      = 0 * K
	lastGraph = Khat
	lastK     = Khat

	mes <- paste( c( iter, " iteration is started.                    " ), collapse = "" )
	cat( mes, "\r" )

	if( save.all == TRUE )
	{
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) )
		{
			result = .C( "bdmcmcExact", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p), 
						allGraphs = as.integer(allGraphs), allWeights = as.double(allWeights), Khat = as.double(Khat), 
						sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
						as.integer(b), as.integer(bstar), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "rjmcmcExact", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p), 
						allGraphs = as.integer(allGraphs), allWeights = as.double(allWeights), Khat = as.double(Khat), 
						sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
						as.integer(b), as.integer(bstar), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) )
		{
			result = .C( "bdmcmcCopula", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						allGraphs = as.integer(allGraphs), allWeights = as.double(allWeights), Khat = as.double(Khat), 
						sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
						as.integer(b), as.integer(bstar), as.double(D), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "rjmcmcCopula", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						allGraphs = as.integer(allGraphs), allWeights = as.double(allWeights), Khat = as.double(Khat), 
						sampleGraphs = as.character(sampleGraphs), graphWeights = as.double(graphWeights), sizeSampleG = as.integer(sizeSampleG),
						as.integer(b), as.integer(bstar), as.double(D), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}	
	}
	else
	{
		if( ( method == "ggm" ) && ( algorithm == "bdmcmc" ) )
		{
			result = .C( "bdmcmcExactPhat", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p), 
						Khat = as.double(Khat), phat = as.double(phat),
						as.integer(b), as.integer(bstar), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}

		if( ( method == "ggm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "rjmcmcExactPhat", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p), 
						Khat = as.double(Khat), phat = as.integer(phat),
						as.integer(b), as.integer(bstar), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}
		
		if( ( method == "gcgm" ) && ( algorithm == "bdmcmc" ) )
		{
			result = .C( "bdmcmcCopulaPhat", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						Khat = as.double(Khat), phat = as.double(phat),
						as.integer(b), as.integer(bstar), as.double(D), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}

		if( ( method == "gcgm" ) && ( algorithm == "rjmcmc" ) )
		{
			result = .C( "rjmcmcCopulaPhat", as.integer(iter), as.integer(burnin), G = as.integer(G), as.double(Ts), K = as.double(K), as.integer(p),
						as.double(Z), as.integer(R), as.integer(n), as.integer(gcgm_NA),
						Khat = as.double(Khat), phat = as.integer(phat),
						as.integer(b), as.integer(bstar), as.double(D), as.double(Ds), as.double(threshold), PACKAGE = "BDgraph" )
		}	
	}

	print( Sys.time() - startTime )  	

	Khat      = matrix( result $ Khat, p, p ) 
	lastGraph = matrix( result $ G, p, p )
	lastK     = matrix( result $ K, p, p )

	colnames( lastGraph ) = colnames( data )

	if( save.all == TRUE )
	{
		if( algorithm == "rjmcmc" ) Khat = Khat / ( iter - burnin )		
		sizeSampleG  = result $ sizeSampleG
		sampleGraphs = result $ sampleGraphs[ 1 : sizeSampleG ]
		graphWeights = result $ graphWeights[ 1 : sizeSampleG ]
		allGraphs    = result $ allGraphs + 1
		allWeights   = result $ allWeights	

		output = list( sampleGraphs = sampleGraphs, graphWeights = graphWeights, Khat = Khat, 
					allGraphs = allGraphs, allWeights = allWeights, lastGraph = lastGraph, lastK = lastK )
	}
	else
	{
		phat   = matrix( result $ phat, p, p ) 
		if( algorithm == "rjmcmc" )
		{
			phat = phat / ( iter - burnin )
			Khat = Khat / ( iter - burnin )
		}
		phat[ lower.tri( phat ) ] = 0
		colnames( phat ) = colnames( data )
		output = list( phat = phat , Khat = Khat, lastGraph = lastGraph, lastK = lastK )
	}
	
	class( output ) = "bdgraph"
	return( output )   
}
      
# summary of bdgraph output
summary.bdgraph = function( object, vis = TRUE, ... )
{
	phat       = object $ phat
	p          = nrow( object $ lastGraph )
	dimlab     = colnames( object $ lastGraph )
	selected_G = matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	

	if( is.null( phat ) )
	{
		sampleGraphs = object $ sampleGraphs
		graphWeights = object $ graphWeights
		max_gWeights = max( graphWeights )
		sum_gWeights = sum( graphWeights )
		max_prob_G   = max_gWeights / sum_gWeights

		if ( is.null( dimlab ) ) dimlab <- as.character( 1 : p )
		vec_G    <- c( rep( 0, p * ( p - 1 ) / 2 ) )		
		indG_max <- sampleGraphs[ which( graphWeights == max_gWeights ) ]
		vec_G[ which( unlist( strsplit( as.character( indG_max ), "" ) ) == 1 ) ] = 1
		selected_G[ upper.tri( selected_G ) ] <- vec_G 
	}
	else
	{
		selected_G[ phat > 0.5 ]  = 1
		selected_G[ phat <= 0.5 ] = 0
	}

	if( vis )
	{
		# plot selected graph (graph with the highest posterior probability)
		G  <- graph.adjacency( selected_G, mode = "undirected", diag = FALSE )
		 
		if( is.null( phat ) ) 
		{
			op = par( mfrow = c( 2, 2 ), pty = "s", omi = c( 0.3, 0.3, 0.3, 0.3 ), mai = c( 0.3, 0.3, 0.3, 0.3 ) ) 
			subGraph = paste( c( "Posterior probability = ", max_prob_G ), collapse = "" )
		}
		else
		{
			subGraph = "Selected graph with edge posterior probability = 0.5"
		}
			
		if( p < 20 ) size = 15 else size = 2
		plot.igraph( G, layout = layout.circle, main = "Selected graph", sub = subGraph, vertex.color = "white", vertex.size = size, vertex.label.color = 'black' )
		
		if( is.null( phat ) )
		{
			# plot posterior distribution of graph
			plot( x = 1 : length( graphWeights ), y = graphWeights / sum_gWeights, type = "h", main = "Posterior probability of graphs",
				 ylab = "Pr(graph|data)", xlab = "graph" )
			
			abline( h = max_prob_G, col = "red" )
			text( which( max_gWeights == graphWeights )[1], max_prob_G, "Pr(selected graph|data)", col = "gray60", adj = c( 0, +1 ) )
			
			# plot posterior distribution of graph size
			sizeSampleGraphs = sapply( sampleGraphs, function(x) length( which( unlist( strsplit( as.character(x), "" ) ) == 1 ) ) )
			xx       <- unique( sizeSampleGraphs )
			weightsg <- vector()

			for( i in 1 : length(xx) ) weightsg[i] <- sum( graphWeights[ which( sizeSampleGraphs == xx[i] ) ] )

			plot( x = xx, y = weightsg / sum_gWeights, type = "h", main = "Posterior probability of graphs size", ylab = "Pr(graph size|data)", xlab = "Graph size" )

			# plot trace of graph size
			allGraphs        = object $ allGraphs
			sizeAllGraphs    = sizeSampleGraphs[ allGraphs ]
			  
			plot( x = 1 : length( allGraphs ), sizeAllGraphs, type = "l", main = "Trace of graph size", ylab = "Graph size", xlab = "Iteration" )
			
			abline( h = sum( selected_G ), col = "red" )	  
			
			par(op)
		}
	}
	
	# phat
	if( is.null( phat ) )
	{
		pvec <- 0 * vec_G
		for( i in 1 : length( sampleGraphs ) )
		{
			which_edge       <- which( unlist( strsplit( as.character( sampleGraphs[i] ), "" ) ) == 1 )
			pvec[which_edge] <- pvec[which_edge] + graphWeights[i]
		}
		phat                  <- 0 * selected_G
		phat[upper.tri(phat)] <- pvec / sum_gWeights
	}
					  
	return( list( selected_G = Matrix( selected_G, sparse = TRUE ), phat = Matrix( phat, sparse = TRUE ), Khat = object $ Khat ) )
}  
   
# plot for class bdgraph
plot.bdgraph = function( x, cut = 0.5, number.g = 1, layout = layout.circle, ... )
{
	phat = x $ phat
	
	if( is.null( phat ) )
	{
		sampleGraphs = x $ sampleGraphs
		graphWeights = x $ graphWeights
		prob_G       = graphWeights / sum( graphWeights )
		sort_prob_G  = sort( prob_G, decreasing = TRUE )

		p            = nrow( x $ lastGraph )
		dimlab       = colnames( x $ lastGraph )
		list_G       = replicate( number.g, matrix( 0, p, p, dimnames = list( dimlab, dimlab ) ), simplify = FALSE )
		vec_G        = c( rep( 0, p * ( p - 1 ) / 2 ) )

		if( number.g == 2 ) op <- par( mfrow = c( 1, 2 ), pty = "s" )
		if( number.g > 2 & number.g < 7 )  op <- par( mfrow = c( 2, number.g %% 2 + trunc( number.g / 2 ) ), pty = "s" )

		for( i in 1 : number.g )
		{
			if( number.g > 6 ) dev.new()  
			indG_i <- sampleGraphs[ which( prob_G == sort_prob_G[i] )[1] ]
			vec_G  <- 0 * vec_G
			vec_G[ which( unlist( strsplit( as.character(indG_i), "" ) ) == 1 ) ] <- 1
			list_G[[i]][ upper.tri( list_G[[i]] ) ] <- vec_G
			G    <- graph.adjacency( list_G[[i]], mode = "undirected", diag = FALSE )
			main <- ifelse( i == 1, "Graph with highest probability", paste( c( i, "th graph" ), collapse = "" ) )
			plot.igraph( G, layout = layout, main = main, sub = paste( c( "Posterior probability = ", 
						round( sort_prob_G[i], 6 ) ), collapse = "" ), ... )	   
		}
		
		if( number.g > 1 & number.g < 7 ) par( op )
	}
	else
	{
		if( ( cut < 0 ) || ( cut > 1 ) ) stop( "Value of 'cut' should be between zero and one." )
		selected_G                = 0 * phat
		selected_G[ phat > cut ]  = 1
		selected_G[ phat <= cut ] = 0		

		G    = graph.adjacency( selected_G, mode = "undirected", diag = FALSE )
		plot.igraph( G, layout = layout, main = "Selected graph", sub = "Edge posterior probability = 0.5", ... )	   		
	}
}
     
# print of the bdgraph output
print.bdgraph = function( x, round = 3, ... )
{
	p_hat = x $ phat
	
	if( is.null( p_hat ) )
	{
		p             = nrow( x $ lastGraph )
		sampleGraphs  = x $ sampleGraphs
		graphWeights  = x $ graphWeights
		# selected graph
		max_gWeights  = max( graphWeights )
		sum_gWeights  = sum( graphWeights )
		vec_G         = c( rep( 0, p * ( p - 1 ) / 2 ) )
		indG_max      = sampleGraphs[ which( graphWeights == max_gWeights )[1] ]
		vec_G[ which( unlist( strsplit( as.character( indG_max ), "" ) ) == 1 ) ] = 1

		dimlab     = colnames( x $ lastGraph )
		selected_G = matrix( 0, p, p, dimnames = list( dimlab, dimlab ) )	
		selected_G[upper.tri(selected_G)] = vec_G
	}
	else
	{
		selected_G                 = 0 * p_hat
		selected_G[ p_hat > 0.5 ]  = 1
		selected_G[ p_hat <= 0.5 ] = 0	
	}
	
	cat( paste( "" ), fill = TRUE )
	cat( paste( "Adjacency matrix of selected graph" ), fill = TRUE )
	cat( paste( "" ), fill = TRUE )
	printSpMatrix( Matrix( selected_G, sparse = TRUE ), col.names = TRUE, note.dropping.colnames = FALSE )

	cat( paste( "" ), fill = TRUE )
	cat( paste( "Size of selected graph = ", sum( selected_G ) ), fill = TRUE )
	if( is.null( p_hat ) )
		cat( paste( "Posterior probability of selected graph = ", max_gWeights / sum_gWeights ), fill = TRUE )  
	else
		cat( paste( "Edge posterior probability of selected graph = ", 0.5 ), fill = TRUE )
	
	cat( paste( "" ), fill = TRUE )
} 
   




