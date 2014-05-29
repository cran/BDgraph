# computing birth and death rates
log_H_ij = function(K, G, i, j, bstar, Ds, p)
{
	e <- c(i, j)
	
	if ( G[i,j] == 0 )
	{
		K12    <- K[e, -e]  
		F      <- K12 %*% solve(K[-e, -e]) %*% t(K12) 

		a      <- K[e, e] - F

		sig    <- sqrt(a[1, 1] / Ds[j, j])
		mu     <- - (Ds[i, j] * a[1, 1]) / Ds[j,j]
		u      <- rnorm(1, mu, sig)
		v      <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[j,j])
		  
		K[i,j] <- u + F[1,2]
		K[j,i] <- u + F[1,2]
		K[j,j] <- v + u ^ 2 / a[1,1] + F[2,2]
	}

	# (i,j) = 0
	K0       <- K
	K0[i, j] <- 0
	K0[j, i] <- 0     
	K_12     <- K0[j, -j, drop = FALSE]
	K0_ij    <- diag(c(K[i, i], K_12 %*% solve(K0[-j, -j]) %*% t(K_12))) 

	# (i,j) = 1
	K_12  <- K[e, -e]  
	K1_ij <- K_12 %*% solve(K[-e, -e]) %*% t(K_12) 

	a11    <- K[i, i] - K1_ij[1, 1]

	log_Hij <- ( (1 / 2) * ( log( Ds[j, j] ) - log( a11 ) )
		        + ((Ds[i, i] - Ds[i, j] ^ 2 / Ds[j, j]) * a11) / 2
		        - sum( Ds[e, e] * (K0_ij - K1_ij) ) / 2  
		       ) # / sqrt(2 * pi)

	return( log_Hij )
}
########################################################################
# sampling from G-Wishart based on Alex Lenkoski (2013) page 122
rgwish.exact <- function( G, b, T, p, threshold = 1e-8 )
{
    nu  <- 1:p
    Psi <- diag( sqrt( rchisq( p, b + p - nu ) ) )

    Psi[upper.tri(Psi)] <- rnorm( p * (p - 1) / 2 )
      
    Psi <- Psi %*% T   
    K   <- t(Psi) %*% Psi 
    
    Sigma      <- solve(K)
    W          <- Sigma
    difference <- 100
    
    while( difference > threshold )
    {
        W.last <- W
        
        for ( j in 1:p )
        {
			# adjacency matrix has to have zero in its diagonal
            N.j  <- (1:p)[G[j, ] == 1]
            
            if ( length(N.j) > 0 )
            {
                W.N.j           <- W[N.j, N.j]
                beta.star.hat.j <- solve(W.N.j, Sigma[N.j, j])
                beta.star       <- rep(0, p)
                beta.star[N.j]  <- beta.star.hat.j             
                ww              <- W[-j, -j] %*% beta.star[-j]
                W[j, -j]        <- ww
                W[-j, j]        <- ww
            } else {
                W[-j, j] <- 0
                W[j, -j] <- 0
            }
        }
          
        difference <- max( abs( W.last - W ) )
    }
      
    K <- solve(W)
    
    return( K )
} 
## Main function: BDMCMC algorithm for selecting the best graphs 
#  based on exchange algorithm to avoid computing the ratio of normalizing constants
#  and exact sampling from precision matrix
bdgraph_ex = function( data, n = NULL, npn = "normal", iter = 5000, 
                       burnin = floor(iter / 2), b = 3, D = NULL, g.start = "full", 
                       K.start = NULL, save.all = FALSE, trace = TRUE )
{
	start.time <- Sys.time()
	
	if ( class(data) == "simulate" ) data <- data $ data

	if ( is.matrix(data) == FALSE & is.data.frame(data) == FALSE ) stop( "Data should be a matrix or dataframe" )
	if ( is.data.frame(data) ) data <- data.matrix(data)
	if ( any( is.na(data) ) ) stop( "Data should contain no missing data" ) 
	if ( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )

	if ( npn != "normal" ) data <- bdgraph.npn(data = data, npn = npn, npn.thresh = NULL)

	dimd <- dim(data)

	if ( isSymmetric(data) )
	{
		if ( is.null(n) ) stop( "Please specify the number of observations 'n'" )
		if ( trace ) cat( "The input is identified as the covriance matrix. \n" )
		S <- data
	} else {
		n <- dimd[1]
		S <- n * cov( data )
		#data <- scale( data )
		#S    <- cor( data )
	}

	p <- dimd[2]

	#if ( g.prior == "Poisson" & is.null(lambda) ) stop( "You should determine value of lambda as the rate of Poisson" )

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

	if ( class(g.start) == "bdgraph" ) 
	{
		G <- g.start $ last.G
		K <- g.start $ last.K
	} else {
		if ( !is.matrix(g.start) )
		{
			if ( g.start == "full" )
			{
				G               <- 0 * S
				G[upper.tri(G)] <- 1
				K               <- matrix( rWishart(1, df = bstar + (p - 2), Sigma = invDs), p, p )
			}

			if ( g.start == "empty" )
			{
				G <- 0 * S
				K <- rgwish.exact(G = G, b = bstar, T = Ts, p = p)
			}

			if ( g.start == "glasso" | g.start == "mb" | g.start == "ct" )
			{
				G               <- huge( data, method = g.start )
				G               <- huge.select(G)
				G               <- as.matrix(G $ refit)
				G[lower.tri(G)] <- 0

				K               <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)
			}
		} else {
			G                         <- as.matrix(g.start)
			G[lower.tri(G, diag = T)] <- 0
			if ( !is.null(K.start) ) K <- K.start else K  <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)  
		}
	}

	sample.G   <- all.G       <- vector() # vector of numbers like "10100"
	weights    <- all.weights <- vector() # waiting time for every state
	sumK       <- 0 * K

	rates <- 0 * K

	for ( g in 1 : iter )
	{
		if ( trace == TRUE && g %% 100 == 0 )
		{
			mes <- paste( c( " iteration: ", g, " from ", iter, ". Graph size= ", sum(G) ), collapse = "" )
			cat(mes, "\r")
			flush.console()	
		}
		
		# using exchange algorithm
		K_prop <- rgwish.exact( G = G + t(G), b = b, T = Ti, p = p )
	
		for ( i in 1 : (p - 1) )
		{
			for ( j in (i + 1) : p )
			{
				logHij <- log_H_ij( K = K, G = G, i = i, j = j, bstar = bstar, Ds = Ds, p = p )
				logI_p <- log_H_ij( K = K_prop, G = G, i = i, j = j, bstar = b, Ds = D, p = p )
				
				rates[i,j] <- if( G[i,j] == 1 ) exp( logHij - logI_p ) else exp( logI_p - logHij )
			
				#if ( g.prior == "Poisson" ) 
				#{
				#	rates[i,j] <- ifelse( G[i,j] == 1, (lambda / (sum(G) + 1)) * rates[i,j], (lambda / (sum(G) + 1)) * rates[i,j] )
				#}
			}
		}

		if ( save.all == TRUE )
		{
			indG        <- paste( G[upper.tri(G)], collapse = '' )
			all.G       <- c( all.G, indG )
			all.weights <- c( all.weights, 1 / sum(rates) )
		}

		if ( g > burnin )
		{
			sumK <- sumK + K
			indG <- paste( G[upper.tri(G)], collapse = '' )
			wh   <- which( sample.G == indG )
			if ( length(weights) != 0 & length(wh) != 0 )
			{
				weights[wh] <- weights[wh] + 1 / sum(rates)
			} else {
				sample.G <- c( sample.G, indG )
				weights  <- c( weights, 1 / sum(rates) )
			}
		}

		# check for numerical issues
		rates.na  <- which( is.na(rates) )
		if ( length(rates.na) > 0)   { rates[rates.na]  <- 0 }
		rates.nan <- which( is.nan(rates) )
		if ( length(rates.nan) > 0)  { rates[rates.nan] <- 0 }  
		# To select new graph
		edge    <- which( rates == max(rates) )[1]
		G[edge] <- 1 - G[edge]

		K <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)
	}

	if ( trace == TRUE )
	{
		mes <- paste( c(" ", iter," iteration done.                    " ), collapse = "" )
		cat( mes, "\r" )
		cat( "\n" )
		flush.console()
		print( Sys.time() - start.time )  
	}

	if ( save.all == TRUE )
	{
		outbdgraph <- list( sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), 
					  all.G = all.G, all.weights = all.weights, last.G = G, last.K = K )
	} else {
		outbdgraph <- list( sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), last.G = G, last.K = K )
	}

	class( outbdgraph ) <- "bdgraph"
	return( outbdgraph )   
}
# copula function: Step 1 in the BDMCMC algorithm
copula = function( Z, K, R, n, p )
{
	for ( j in sample(1:p) ) 
	{   # looping randomly through the variables (= columns of Y)
		sdj <- sqrt( 1 / K[j, j] )     # 2a: # variance of component j (given the rest!)
		# Kjc <- -K[j, -j] / K[j, j]     # Hoff page 273 part of 2b
		# muj <- Z[ , -j] %*% t(t(Kjc))  # 2b
		# matrixproduct for (n x p-1) * (p-1 x 1) = (n x 1)
		muj <- Z[ ,-j, drop = FALSE] %*% - K[-j,j, drop = FALSE] / K[j,j]	 # Hoff page 273 part of 2b

		# interval and sampling
		for( r in sort( unique(R[ , j]) ) )
		{
			ir      <- (1:n)[R[ , j] == r & !is.na(R[ , j])]
			lb      <- suppressWarnings( max( Z[ R[ , j] < r, j], na.rm = T) )
			ub      <- suppressWarnings( min( Z[ R[ , j] > r, j], na.rm = T) )
			# truncated normal distr.
			Z[ir,j] <- qnorm( runif( length(ir), pnorm(lb, muj[ir], sdj), pnorm(ub, muj[ir], sdj) ), muj[ir], sdj )
		}

		ir       <- (1:n)[is.na(R[ , j])]
		Z[ir, j] <- rnorm( length(ir), muj[ir], sdj )
	}

	return( Z )
}
# BDMCMC algorithm with copula approach of mixed data
bdgraph_copula = function( data, n = NULL, iter = 5000, burnin = floor(iter / 2), 
                           b = 3, D = NULL, g.start = "full", K.start = NULL, 
				           save.all = FALSE, trace = TRUE )
{
	start.time <- Sys.time()
	if (class(data) == "simulate") data <- data $ data

	if ( is.matrix(data) == FALSE & is.data.frame(data) == FALSE ) stop( "Data should be a matrix or dataframe" )
	if ( is.data.frame(data) ) data <- data.matrix(data)
	#if ( any( is.na(data) ) ) stop( "Data should contain no missing data" ) 
	if ( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )

	#if ( npn != "normal" ) data <- bdgraph.npn(data = data, npn = npn, npn.thresh = NULL)

	dimd <- dim(data)
	p    <- dimd[2]
	n    <- dimd[1]

	# only for copula
	Z              <- qnorm( apply(data, 2, rank, ties.method = "random") / (n + 1) )
	Zfill          <- matrix(rnorm(n * p), n, p)   # for missing values
	Z[is.na(data)] <- Zfill[is.na(data)]        # for missing values
	Z              <- t( (t(Z) - apply(Z, 2, mean)) / apply(Z, 2, sd) )
	S              <- t(Z) %*% Z

#	if ( isSymmetric(data) )
#	{
#		if ( is.null(n) ) stop( "Please specify the number of observations 'n'" )
#		if ( trace ) cat( "The input is identified as the covriance matrix. \n" )
#		S <- data
#	} else {
#		data <- scale( data )
#		S    <- cor( data )
#	}

#	if ( g.prior == "Poisson" & is.null(lambda) ) stop( "You should determine value of lambda as the rate of Poisson" )

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

	if ( class(g.start) == "bdgraph" ) 
	{
		G <- g.start $ last.G
		K <- g.start $ last.K
	} else {
		if ( !is.matrix(g.start) )
		{
			if ( g.start == "full" )
			{
				G               <- 0 * S
				G[upper.tri(G)] <- 1
				K               <- matrix( rWishart(1, df = bstar + (p - 2), Sigma = invDs), p, p )
			}

			if ( g.start == "empty" )
			{
				G <- 0 * S
				K <- rgwish.exact(G = G, b = bstar, T = Ts, p = p)
			}

			if ( g.start == "glasso" | g.start == "mb" | g.start == "ct" )
			{
				G               <- huge( data, method = g.start )
				G               <- huge.select(G)
				G               <- as.matrix(G $ refit)
				G[lower.tri(G)] <- 0

				K               <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)
			}
		} else {
			G                         <- as.matrix(g.start)
			G[lower.tri(G, diag = T)] <- 0
			if ( !is.null(K.start) ) K <- K.start else K  <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)  
		}
	}

	sample.G <- all.G       <- vector() # vector of numbers like "10100"
	weights  <- all.weights <- vector() # waiting time for every state
	sumK     <- 0 * K

	# only for copula
	R <- NULL
	for(j in 1:p) { R <- cbind(R, match(data[ , j], sort(unique(data[ , j])))) }

	rates <- 0 * K

	for ( g in 1 : iter )
	{
		if ( trace == TRUE && g %% 100 == 0 )
		{
			mes <- paste( c( " iteration: ", g, " from ", iter, ". Graph size= ", sum(G) ), collapse = "" )
			cat(mes, "\r")
			flush.console()	
		}
		
		# First step: copula 
		# here we will use a copula function
		Z <- copula(Z = Z, K = K, R = R, n = n, p = p)
		S <- t(Z) %*% Z

		Ds    <- D + S

		# Second step: BD-MCMC iteration  		
		
		# using exchange algorithm
		K_prop <- rgwish.exact( G = G + t(G), b = b, T = Ti, p = p )
	
		for ( i in 1 : (p - 1) )
		{
			for ( j in (i + 1) : p )
			{
				logHij <- log_H_ij( K = K, G = G, i = i, j = j, bstar = bstar, Ds = Ds, p = p )
				logI_p <- log_H_ij( K = K_prop, G = G, i = i, j = j, bstar = b, Ds = D, p = p )
				
				rates[i,j] <- if( G[i,j] == 1 ) exp( logHij - logI_p ) else exp( logI_p - logHij )
			
				#if ( g.prior == "Poisson" ) 
				#{
				#	rates[i,j] <- ifelse( G[i,j] == 1, (lambda / (sum(G) + 1)) * rates[i,j], (lambda / (sum(G) + 1)) * rates[i,j] )
				#}
			}
		}

		if ( save.all == TRUE )
		{
			indG        <- paste( G[upper.tri(G)], collapse = '' )
			all.G       <- c( all.G, indG )
			all.weights <- c( all.weights, 1 / sum(rates) )
		}

		if ( g > burnin )
		{
			sumK <- sumK + K
			indG <- paste( G[upper.tri(G)], collapse = '' )
			wh   <- which( sample.G == indG )
			if ( length(weights) != 0 & length(wh) != 0 )
			{
				weights[wh] <- weights[wh] + 1 / sum(rates)
			} else {
				sample.G <- c( sample.G, indG )
				weights  <- c( weights, 1 / sum(rates) )
			}
		}

		# check for numerical issues
		rates.na  <- which( is.na(rates) )
		if ( length(rates.na) > 0)   { rates[rates.na]  <- 0 }
		rates.nan <- which( is.nan(rates) )
		if ( length(rates.nan) > 0)  { rates[rates.nan] <- 0 }  
		# To select new graph
		edge    <- which( rates == max(rates) )[1]
		G[edge] <- 1 - G[edge]

		K <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)
	}

	if ( trace == TRUE )
	{
		mes <- paste( c(" ", iter," iteration done.                    " ), collapse = "" )
		cat( mes, "\r" )
		cat( "\n" )
		flush.console()
		print( Sys.time() - start.time )  
	}

	if ( save.all == TRUE )
	{
		outbdgraph <- list( sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), 
					  all.G = all.G, all.weights = all.weights, last.G = G, last.K = K )
	} else {
		outbdgraph <- list( sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), last.G = G, last.K = K )
	}

	class( outbdgraph ) <- "bdgraph"
	return( outbdgraph )   
}
# computing birth and death rates
rate_ij = function(K, A, i, j, b, bstar, Ds, Ti, p)
{
	e          <- c(i, j)
	current_ij <- A[i,j]
	
	if (A[i,j] == 0)
	{
		A[i,j] <- 1

		K12    <- K[e, -e]  
		F      <- K12 %*% solve(K[-e, -e]) %*% t(K12) 

		a      <- K[e, e] - F

		sig    <- sqrt(a[1, 1] / Ds[j, j])
		mu     <- - (Ds[i, j] * a[1, 1]) / Ds[j,j]
		u      <- rnorm(1, mu, sig)
		v      <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[j,j])
		  
		K[i,j] <- u + F[1,2]
		K[j,i] <- u + F[1,2]
		K[j,j] <- v + u ^ 2 / a[1,1] + F[2,2]
	}

	# (i,j) = 0
	K0       <- K
	K0[i, j] <- 0
	K0[j, i] <- 0     
	K_12     <- K0[j, -j, drop = FALSE]
	K0_ij    <- diag(c(K[i, i], K_12 %*% solve(K0[-j, -j]) %*% t(K_12))) 

	# (i,j) = 1
	K_12  <- K[e, -e]  
	K1_ij <- K_12 %*% solve(K[-e, -e]) %*% t(K_12) 

	a11    <- K[i, i] - K1_ij[1, 1]
	nustar <- sum(A[i, ])

	rate <- sqrt(2) * Ti[i, i] * Ti[j, j] *
		    exp( lgamma((b + nustar) / 2) - lgamma((b + nustar - 1) / 2) 
		    + (1 / 2) * ( log( Ds[j, j] ) - log( a11 ) )
		    + ((Ds[i, i] - Ds[i, j] ^ 2 / Ds[j, j]) * a11) / 2
		    - sum( Ds[e, e] * (K0_ij - K1_ij) ) / 2 )

	if ( current_ij == 0 ) rate <- 1 / rate 

	if ( !is.finite(rate) ) rate <- gamma(171)

	return(rate)
}
# Main function: BDMCMC algorithm for graph determination
# based on our approximation for computing the ratio of normalizing constants
#  and exact sampling from precision matrix
bdgraph_ap = function( data, n = NULL, npn = "normal", iter = 5000, 
                       burnin = floor(iter / 2), b = 3, D = NULL, g.start = "full", 
                       K.start = NULL, save.all = FALSE, trace = TRUE )
{
	start.time <- Sys.time()
	if (class(data) == "simulate") data <- data $ data

	if ( is.matrix(data) == FALSE & is.data.frame(data) == FALSE ) stop( "Data should be a matrix or dataframe" )
	if ( is.data.frame(data) ) data <- data.matrix(data)
	if ( any(is.na(data)) ) stop( "Data should contain no missing data" ) 
	if ( iter <= burnin )   stop( "Number of iteration must be more than number of burn-in" )

	if ( npn != "normal" ) data <- bdgraph.npn( data = data, npn = npn, npn.thresh = NULL )

	dimd <- dim(data)

	if ( isSymmetric(data) )
	{
		if ( is.null(n) ) stop( "Please specify the number of observations 'n'" )
		if ( trace ) cat( "The input is identified as the covriance matrix. \n" )
		S <- data
	} else {
		n <- dimd[1]
		S <- n * cov( data )
	}

	p <- dimd[2]

	#if ( g.prior == "Poisson" & is.null(lambda) ) stop( "You should determine value of lambda as the rate of Poisson" )

	if ( is.null(D) )
	{ 
		D  <- diag(p)
		Ti <- D
		H  <- D
	} else {
		Ti <- chol( solve(D) )
		H  <- Ti / t( matrix( rep(diag(Ti) ,p), p, p ) )
	}

	bstar <- b + n
	Ds    <- D + S
	invDs <- solve(Ds)
	Ts    <- chol(invDs)
	Hs    <- Ts / t( matrix( rep(diag(Ts) ,p), p, p ) )	

	if ( class(g.start) == "bdgraph" ) 
	{
		G <- g.start $ last.G
		K <- g.start $ last.K
	} else {
		if ( !is.matrix(g.start) )
		{
			if ( g.start == "full" )
			{
				G               <- 0 * S
				G[upper.tri(G)] <- 1
				K               <- matrix(rWishart(1, df = bstar + (p - 2), Sigma = invDs), p, p)
			}

			if ( g.start == "empty" )
			{
				G <- 0 * S
				K <- rgwish.exact(G = G, b = bstar, T = Ts, p = p)
			}

			if ( g.start == "glasso" | g.start == "mb" | g.start == "ct" )
			{
				G               <- huge(data, method = g.start)
				G               <- huge.select(G)
				G               <- as.matrix(G $ refit)
				G[lower.tri(G)] <- 0

				K               <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)
			}
		} else {
			G                          <- as.matrix(g.start)
			G[lower.tri(G, diag = T)]  <- 0
			if ( !is.null(K.start) ) K <- K.start else K  <- rgwish.exact(G = G + t(G), b = bstar, T = Ts, p = p)  
		}
	}

	sample.G   <- all.G       <- vector() # vector of numbers like "10100"
	weights    <- all.weights <- vector() # waiting times for each state
	sumK       <- 0 * K

	rates <- 0 * K

	for ( g in 1 : iter )
	{
		if ( trace == TRUE && g %% 100 == 0 )
		{
			mes <- paste( c( " iteration: ", g, " from ", iter, ". Graph size= ", sum(G) ), collapse = "" )
			cat(mes, "\r")
			flush.console()	
		}
	
		for ( i in 1 : (p - 1) )
		{
			for ( j in (i + 1) : p )
			{
				rates[i,j] <- rate_ij(K = K, A = G, i = i, j = j, b = b, bstar = bstar, Ds = Ds, Ti = Ti, p = p)
			
				#if ( g.prior == "Poisson" ) 
				#{
				#	rates[i,j] <- ifelse( G[i,j] == 1, (lambda / (sum(G) + 1)) * rates[i,j], (lambda / (sum(G) + 1)) * rates[i,j] )
				#}
			}
		}

		if ( save.all == TRUE )
		{
			indG        <- paste( G[upper.tri(G)], collapse = '' )
			all.G       <- c( all.G, indG )
			all.weights <- c( all.weights, 1 / sum(rates) )
		}

		if ( g > burnin )
		{
			sumK <- sumK + K
			indG <- paste( G[upper.tri(G)], collapse = '' )
			wh   <- which(sample.G == indG)
			if ( length(weights) != 0 & length(wh) != 0 )
			{
				weights[wh] <- weights[wh] + 1 / sum(rates)
			} else {
				sample.G <- c( sample.G, indG )
				weights  <- c( weights, 1 / sum(rates) )
			}
		}
		
		# check for numerical issues
		rates.inf <- which( is.infinite(rates) )
		if ( length(rates.inf) > 0 ) { rates[rates.inf] <- gamma(170); warning( "Nomerical issue: some of rates are infinite" ) }
		rates.na  <- which( is.na(rates) )
		if ( length(rates.na) > 0)   { rates[rates.na]  <- 0; warning( "Nomerical issue: some of rates are NA" ) }
		rates.nan <- which( is.nan(rates) )
		if ( length(rates.nan) > 0)  { rates[rates.nan] <- 0; warning( "Nomerical issue: some of rates are NaN" ) }  
		# To select new graph
		if ( sum(rates) > 0 )
		{
			edge    <- which( rmultinom( 1, 1, rates ) == 1 )
			G[edge] <- 1 - G[edge]
		} else {
			warning( "Nomerical issue: Sum of rates is zero" )
		}
		 
		if ( p < 11 )
		{
			K <- block.gibbs.low(K = K, A = G, bstar = bstar, Ds = Ds, p = p)
		} else { 
			K <- block.gibbs.high(K = K, A = G, bstar = bstar, Ds = Ds, p = p)
		}
	}

	if ( trace == TRUE )
	{
		mes <- paste(c(" ", iter," iteration done.                    "), collapse = "")
		cat( mes, "\r" )
		cat("\n")
		flush.console()
		print( Sys.time() - start.time )  
	}

	if ( save.all == TRUE )
	{
		outbdgraph <- list( sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), 
					  all.G = all.G, all.weights = all.weights, last.G = G, last.K = K )
	} else {
		outbdgraph <- list( sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), last.G = G, last.K = K )
	}

	class( outbdgraph ) <- "bdgraph"
	return( outbdgraph )   
}
## Main function: BDMCMC algorithm for selecting the best graphs 
bdgraph = function( data, n = NULL, method = "exact", npn = "normal", 
                    iter = 5000, burnin = floor(iter / 2), b = 3, D = NULL, 
                    g.start = "full", K.start = NULL, save.all = FALSE, trace = TRUE )
{
	if ( method == "exact" )
	{
		output <- bdgraph_ex( data = data, n = n, npn = npn, 
							  iter = iter, burnin = burnin, b = b, D = D, 
                              g.start = g.start, K.start = K.start, 
			                  save.all = save.all, trace = trace )
	}
	
	if ( method == "approx" )
	{
		output <- bdgraph_ap( data = data, n = n, npn = npn, 
                              iter = iter, burnin = burnin, b = b, D = D, 
                              g.start = g.start, K.start = K.start, 
			                  save.all = save.all, trace = trace )
	}

	if ( method == "copula" )
	{
		output <- bdgraph_copula( data = data, n = n, iter = iter, burnin = burnin, 
                                  b = b, D = D, g.start = g.start, K.start = K.start, 
				                  save.all = save.all, trace = trace )
	}
	
	return( output )
}
# computing probability of all links of the graph
phat = function(output, round = 3)
{
	sample.G <- output $ sample.G
	weights  <- output $ weights
	p        <- nrow(output $ last.G)
	pvec     <- c(rep(0, p * (p - 1) / 2))

	for (i in 1 : length(sample.G))
	{
		inp       <- which(unlist(strsplit(as.character(sample.G[i]), "")) == 1)
		pvec[inp] <- pvec[inp] + weights[i]
	}

	dimlab <- dimnames(output $ last.G)
	if (is.null(dimlab))
	{
		dimlab <- as.character(1 : p)
		phat   <- matrix(0, p, p, dimnames = list(dimlab, dimlab))
	} else {
		phat   <- matrix(0, p, p, dimnames = dimlab)
	}

	phat[upper.tri(phat)] <- pvec / sum(weights)

	return(Matrix(round(phat, round)))
}
# To check the convergency of the BDMCMC algorithm
plotcoda = function(output, thin = NULL, trace = TRUE, main = NULL, ...)
{
	if (is.null(output $ all.G)) stop("Function needs output of 'bdgraph' with option save.all = T")  
	if (is.null(thin)) thin = ceiling(length(output $ all.G) / 1000)

	op          <- par(mfrow = c(2, 2), pty = "s")
	p           <- nrow(output $ last.G)
	all.weights <- output $ all.weights
	all.G       <- output $ all.G

	weights <- output $ weights
	bestg   <- output $ sample.G[which(max(weights) == weights)]	
	lin     <- length(which(unlist(strsplit(as.character(bestg), "")) == 1))
	y       <- sapply(all.G, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))

	plot(x = (1 : length(all.G)), y, type = "l", main = "Trace of graph size",
	   ylab = "graph size", xlab = "iteration", ...)
	abline(h = lin, col = "red")
	acf(y, main = "ACF for graph size")
	pacf(y, main = "PACF for graph size")

	allG.new        <- all.G[c(thin * (1 : floor(length(all.G) / thin)))]
	all.weights.new <- all.weights[c(thin * (1 : floor(length(all.weights) / thin)))]
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
		ffv[inp] <- ffv[inp] + all.weights.new[g]
		ff[ ,g]  <- ffv / sum(all.weights.new[c(1 : g)])    	 
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
    if (is.null(output $ all.G)) stop("This function needs output of 'bdgraph' with option save.all = T")
	
    all.G   <- output $ all.G
	weights <- output $ weights
	bestg   <- output $ sample.G[which(max(weights) == weights)]	
	lin     <- length(which(unlist(strsplit(as.character(bestg), "")) == 1))
    y       <- sapply(all.G, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
	
	if (is.null(main)) main = "Trace of graph size"
	
	if (acf == FALSE & pacf == FALSE)
	{
		plot(x = 1 : length(all.G), y, type = "l", main = main,
			ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	   
	}
	
	if (acf == TRUE & pacf == TRUE)
	{
		op <- par(mfrow = c(2, 2), pty = "s") 
		plot(x = 1 : length(all.G), y, type = "l", main = main,
			   ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	  
		acf(y,  main = "ACF for graph size")
		pacf(y, main = "PACF for graph size")
		par(op)
	}
	
	if (acf == TRUE & pacf == FALSE)
	{
		op <- par(mfrow = c(1, 2), pty = "s") 
		plot(x = 1 : length(all.G), y, type = "l", main = main,
			   ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	  
		acf(y, main = "ACF for graph size")
		par(op)
	}
	
	if (acf == FALSE & pacf == TRUE)
	{
		op <- par(mfrow = c(1, 2), pty = "s") 
		plot(x = 1 : length(all.G), y, type = "l", main = main,
			   ylab = "graph size", xlab = "iteration", ...)
		abline(h = lin, col = "red")	  
		pacf(y, main = "PAIC for graph size")
		par(op)
	}		
}  
# To select the best graph (graph with the highest posterior probability) 
select = function (output, vis = FALSE)
{
	sample.G   <- output $ sample.G
	weights    <- output $ weights
	p          <- nrow(output $ last.G)
	prob.G     <- weights / sum(weights)
	max.prob.G <- which(prob.G == max(prob.G))
	
	if(length(max.prob.G) > 1) max.prob.G <- max.prob.G[1] 
	
	gi        <- sample.G[max.prob.G]
	gv        <- c(rep(0, p * (p - 1) / 2))
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1

	dimlab   <- dimnames(output $ last.G)
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
	sample.G <- output $ sample.G
	weights  <- output $ weights

	if (is.null(G))
	{
		p      <- nrow(output $ last.G)
		graphi <- list()
		gv     <- c(rep(0, p * (p - 1) / 2))  

		for (i in 1 : g)
		{
			gi <- sample.G[which(weights == sort(weights, decreasing = T)[i])]
			gv <- 0 * gv
			gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
			graphi[[i]] <- matrix(0, p, p)
			graphi[[i]][upper.tri(graphi[[i]])] <- gv
			graphi[[i]] <- Matrix(graphi[[i]] + t(graphi[[i]]), sparse = TRUE)
		}

		return(list(best.G = graphi, prob.G = sort(weights, decreasing = T)[1 : g] / sum(weights)))
		
	} else {
		if (class(G) == "simulate") G <- G $ G

		G     <- as.matrix(G)
		indA  <- paste(G[upper.tri(G)], collapse = '')
		wh    <- which(sample.G == indA)
		probG <- ifelse(length(wh) == 0, 0, weights[wh] / sum(weights))

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
                        K = NULL, sigma = NULL, v = 0.3, u = 0.2, mean = 0, 
                        vis = FALSE )
{
    if ( is.matrix(K) | is.matrix(K) )  graph <- "fixed"
    
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
	    G       <- toeplitz(c(0, 1, rep(0, p - 2)))
        G[1, p] <- 1
		G[p, 1] <- 1
	}

	if ( graph == "scale-free" )
	{
		data_sim <- huge.generator( n = 2, d = 5, graph = "scale-free" )
		G        <- as.matrix( data_sim $ theta ) 
	}
	
    if ( !is.null(sigma) ) K <- solve(sigma)   

    if ( is.matrix(K) )
    { 
		G     <- 1 * ( abs(K) > 0.02 )
		if( is.null(sigma) ) sigma <- solve(K)
		
    } else {
		if ( !is.null(K) && K == "gwish" )
		{
			Ti      <- chol( solve(D) )
			diag(G) <- 0
			K       <- rgwish.exact( G = G, b = b, T = Ti, p = p )
			sigma   <- solve(K)
		}
	}

    if ( is.null(sigma) & is.null(K) )
	{    	
		K       <- G * v
		# make K positive definite and standardized
		diag(K) <- abs(min(eigen(K) $ values)) + u
		sigma   <- cov2cor(solve(K))
		K       <- solve(sigma)
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
	list.G  <- x $ sample.G
	weights <- x $ weights
	p       <- nrow(x $ last.G)
	prob.G  <- weights / sum(weights)
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
		dimnames(graphi[[i]]) <- dimnames(x $ last.G)
		G    <- graph.adjacency(graphi[[i]], mode = "undirected", diag = FALSE)
		main <- ifelse (i == 1, "Graph with highest probability", paste(c(i, "th graph"), collapse = ""))
		plot.igraph(G, layout = layout, main = main, sub = paste(c("Posterior probability = ", 
		            round(sort(prob.G, decreasing = TRUE)[i], 4)), collapse = ""), ...)	   
	}
	
	if (g > 1 & g < 7) par(op)
}
# summary of bdgraph output
summary.bdgraph = function(object, vis = TRUE, layout = layout.circle, ...)
{
	sample.G <- object $ sample.G
	weights  <- object $ weights
	p        <- nrow(object $ last.G)
	gv       <- c(rep(0, p * (p - 1) / 2))

	dimlab   <- dimnames(object $ last.G)
	if (is.null(dimlab))
	{ 
		dimlab <- as.character(1 : p)
		graphi <- matrix(0, p, p, dimnames = list(dimlab, dimlab))
	} else {
		graphi <- matrix(0, p, p, dimnames = dimlab)
	}	

	prob.G <- weights / sum(weights)
	gi     <- sample.G[which(prob.G == max(prob.G))]
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
	graphi[upper.tri(graphi)] <- gv 

	if (vis == TRUE)
	{
		# plot best graph
		G  <- graph.adjacency(graphi, mode = "undirected", diag = FALSE)
		 
		op <- par(mfrow = c(2, 2), pty = "s")

		plot.igraph(G, layout = layout, main = "Best graph",
		  sub = paste(c("Posterior probability = ", round(max(prob.G), 4)), collapse = ""), ...)
		# plot posterior distribution of graph
		plot(x = 1 : length(weights), y = weights / sum(weights), type = "h", main = "Posterior probability",
			 ylab = "Pr(graph|data)", xlab = "graph")
		abline(h = max(weights) / sum(weights), col = "red")
		text(which(max(weights) == weights), max(weights) / sum(weights), "P(best graph|data)", col = "gray60", adj = c(0, + 1))
		# plot posterior distribution of graph size
		suma     <- sapply(sample.G, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
		xx       <- unique(suma)
		weightsg <- vector()

		for (i in 1 : length(xx))
		{
			weightsg[i] <- sum(weights[which(suma == xx[i])])
		}

		plot(x = xx, y = weightsg / sum(weights), type = "h", main = "Posterior probability",
			 ylab = "Pr(graph size|data)", xlab = "graph size")

		if (!is.null(object $ all.G))
		{
			# plot trace of graph size
			if (!is.null(object $ all.G))
			{
				all.G <- object $ all.G
				yy    <- sapply(all.G, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
				plot(x = 1 : length(all.G), yy, type = "l", main = "Trace for graph size",
				  ylab = "graph size", xlab = "iteration")
				abline(h = sum(graphi), col = "red")	  
			}
		}
		par(op)
	}
	# phat
	pvec <- 0 * gv
	for (i in 1 : length(sample.G))
	{
		inp       <- which(unlist(strsplit(as.character(sample.G[i]), "")) == 1)
		pvec[inp] <- pvec[inp] + weights[i]
	}

	phat                  <- 0 * graphi
	phat[upper.tri(phat)] <- pvec / sum(weights)
	# estimation for precision matrix 
	Khat        <- object $ Khat

	return.list <- list(best.graph = Matrix(graphi + t(graphi), sparse = TRUE), phat = Matrix(round(phat, 2), 
					  sparse = TRUE), Khat = round(Khat, 3))
					  
	return(return.list)
}
# print of the bdgraph output
print.bdgraph = function(x, round = 3, Khat = FALSE, phat = FALSE, ...)
{
	sample.G <- x $ sample.G
	weights  <- x $ weights
	p        <- nrow(x $ last.G)
	# best graph
	prob.G   <- weights / sum(weights)
	gv       <- c(rep(0, p * (p - 1) / 2))
	gi       <- sample.G[which(prob.G == max(prob.G))]
	gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1

	dimlab   <- dimnames(x $ last.G)
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
	cat(paste("Posterior probability of best graph = ", round(max(weights) / sum(weights), round)), fill = TRUE)  
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
		for (i in 1 : length(sample.G))
		{
			inp       <- which(unlist(strsplit(as.character(sample.G[i]), "")) == 1)
			pvec[inp] <- pvec[inp] + weights[i]
		}

		phat                  <- 0 * graphi
		phat[upper.tri(phat)] <- pvec / sum(weights)
		cat(paste(""), fill = TRUE)
		cat(paste("Posterior probability of links"), fill = TRUE)
		cat(paste(""), fill = TRUE)

		printSpMatrix(Matrix(round(phat, round), sparse = TRUE), col.names = TRUE, note.dropping.colnames = FALSE)  
	}
} 
# sampling from G-Wishart distribution
rgwish = function( n = 1, G = NULL, b = 3, D = NULL, method = "exact", start.K = NULL )
{
	if (is.null(G)) stop("You should determine the adjacency matrix G")
	G <- as.matrix(G)
	if (sum((G == 1) * (G == 0)) != 0) stop("Elements of matrix G should be zero or one")

	if (sum(upper.tri(G)) == sum(G[upper.tri(G == 1)])) method <- "rWishart"
	
	G[lower.tri(G, diag(TRUE))] <- 0
	p <- nrow(G)  
	
	if ( is.null(D) ) D <- diag(p)	

	samples <- array(0, c(p, p, n))

	if (method == "exact")
	{
		for (i in 1 : n)
		{
			Ti           <- chol( solve(D) ) 
			samples[,,i] <- rgwish.exact(G = G + t(G), b = b, T = Ti, p = p, threshold = 1e-8)
		}	
	}

	if ( method == "block gibbs" ) 
	{
		if ( is.null(start.K) )
		{
			start.K <- diag(p)
		}

		for (i in 1 : n)
		{
			if (p < 11)
			{
				samples[,,i] <- block.gibbs.low (K = start.K, A = G, bstar = b, Ds = D, p = p)
			} else { 
				samples[,,i] <- block.gibbs.high (K = start.K, A = G, bstar = b, Ds = D, p = p)
			}

			start.K      <- samples[,,i]
		}
	}

	if ( method == "accept-reject") 
	{
		Ts <- chol(solve(D))
		H  <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))

		for (i in 1 : n)
		{
			samples[,,i] <- sampleK(A = G, b = b, H = H, Ts = Ts, p = p, iter = 1)
		}
	}

	if ( method == "rWishart" ) samples <- rWishart(n = n, df = b + p - 1, Sigma = solve(D))

	return(samples)   
}
# Accept-reject algorithm: sampling from precision matrix
sampleK = function(A, b, H, Ts, p, iter = 1)
{
	psi  <- Psi( A, b, H, p )
	cont <- 0

	for ( i in 1 : (iter * 1000) )
	{  
		psi.new <- Psi(A, b, H, p)
		alpha   <- exp((sum((1 - (diag(p) + A)) * psi * psi) - 
		           sum((1 - (diag(p) + A)) * psi.new * psi.new)) / 2)

		if (is.nan(alpha)) alpha <- 0

		if (runif(1) < alpha)
		{ 
			cont    <- cont + 1
			psi.new <- psi
		} 

		if (cont == iter) break
	}

	sai <- psi %*% Ts
	return( t(sai) %*% sai )
}
# Blocked Gibbs algorithm for sampling from G-Wishart distribution 
# pairwise blocked Gibbs algorithm for low-dimensional graphs 
block.gibbs.low = function (K, A, bstar, Ds, p)
{
	sumcol <- colSums(A)
	sumrow <- rowSums(A)

	for (i in 1 : (p - 1))
	{
		if (sumrow[i] != 0)
		{
			for (j in (i + 1) : p)
			{
				if (A[i, j] == 1)
				{
					pair <- c(i, j)
					B    <- Ds[pair, pair]
					a    <- rWishart(1, df = bstar + 1, Sigma = solve(B))
					k12  <- K[pair, - pair]
					Kc   <- matrix(a, 2, 2) + k12 %*% solve(K[- pair, - pair]) %*% t(k12)
					K[pair, pair] <- (Kc + t(Kc)) / 2
				}
			}
		}

		if (sumrow[i] + sumcol[i] == 0)
		{
			k12     <- K[i, - i, drop = FALSE]
			K[i, i] <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[i, i]) + 
					 k12 %*% solve(K[- i, - i, drop = FALSE]) %*% t(k12)
		}
	}

	if (sumcol[p] == 0)
	{
		k12     <- K[p, - p, drop = FALSE]
		K[p, p] <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[p, p]) + 
			   k12 %*% solve(K[- p, - p, drop = FALSE]) %*% t(k12)
	}

	return(K)
}
# pairwise blocked Gibbs sampler algorithm for high-dimensional graphs
block.gibbs.high = function (K, A, bstar, Ds, p)
{
	sumcol <- colSums(A)
	sumrow <- rowSums(A)
	Sig    <- solve(K)

	for (i in 1 : (p - 1))
	{
		if (sumrow[i] != 0)
		{
			for (j in (i + 1) : p)
			{
				if (A[i, j] == 1)
				{
					pair     <- c(i, j)
					B        <- Ds[pair, pair]
					a        <- rWishart(1, df = bstar + 1, Sigma = solve(B))
					k12      <- K[pair, - pair]
					Sig12    <- Sig[pair, - pair]
					Sig22    <- Sig[- pair, - pair] 
					invSig11 <- solve(Sig[pair, pair])
					invSig11 <- (invSig11 + t(invSig11)) / 2 # Numerical stable
					invk22   <- Sig22 - t(Sig12) %*% invSig11 %*% Sig12
					Kc       <- matrix(a, 2, 2) + k12 %*% invk22 %*% t(k12)
					Kc       <- (Kc + t(Kc)) / 2
					Delta    <- solve(K[pair, pair] - Kc)
					K[pair, pair] <- Kc
					# step 2: update Sigma 
					Sigbb    <- Sig[pair, pair]
					aa       <- solve(Delta - Sigbb)
					aa       <- (aa + t(aa)) / 2
					Sig      <- Sig + Sig[ , pair] %*% aa %*% t(Sig[ , pair])
				}
			}
		}

		if (sumrow[i] + sumcol[i] == 0)
		{
			a        <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[i, i])
			k12      <- K[i, - i, drop = FALSE]
			Sig12    <- Sig[i, - i, drop = FALSE]
			Sig22    <- Sig[- i, - i, drop = FALSE] 
			invSig11 <- solve(Sig[i, i])
			invSig11 <- (invSig11 + t(invSig11)) / 2 # Numerical stable
			invk22   <- Sig22 - t(Sig12) %*% invSig11 %*% Sig12
			Kc       <- a + k12 %*% invk22 %*% t(k12)
			Delta    <- solve(K[i, i] - Kc)
			K[i, i]  <- Kc
			# step 2: update Sigma 
			Sigbb    <- Sig[i, i]
			aa       <- solve(Delta - Sigbb)
			Sig      <- Sig + Sig[ , i, drop = FALSE] %*% aa %*% t(Sig[ , i, drop = FALSE])	  
		}
	}

	if (sumcol[p] == 0)
	{
		a        <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[p, p])
		k12      <- K[p, - p, drop = FALSE]
		Sig12    <- Sig[p, - p, drop = FALSE]
		Sig22    <- Sig[- p, - p, drop = FALSE] 
		invSig11 <- solve(Sig[p, p])
		invSig11 <- (invSig11 + t(invSig11)) / 2 # Numerical stable
		invk22   <- Sig22 - t(Sig12) %*% invSig11 %*% Sig12
		Kc       <- a + k12 %*% invk22 %*% t(k12)
		Delta    <- solve(K[p, p] - Kc)
		K[p, p]  <- Kc
		# step 2: update Sigma 
		Sigbb    <- Sig[p, p]
		aa       <- solve(Delta - Sigbb)
		Sig      <- Sig + Sig[ , p, drop = FALSE] %*% aa %*% t(Sig[ , p, drop = FALSE])		
	}

	return(K)
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
# To simulate the Psi matrix 
Psi = function(A, b, H, p)
{
	nu          <- rowSums(A)
	psi         <- diag(sqrt(rchisq(p, b + nu)))
	psi[A == 1] <- rnorm(1)

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
			}
		}
	}

	return(psi)
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






