## All functions for the "BDgraph" package
###########################################################
## function for only obtaining elements of Psi matrix
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
# Algorithm 3.1: for Monte Carlo approximation for expectation in normalizing constant
Exp.MC = function(A, b, H, mc, p)
{
   nu  <- rowSums(A)  
   f_T <- c(rep(0, mc))
   
   for (k in 1 : mc)
   {
      # calculating ratio of expectation for birth/death rates 
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
	
    return( mean( exp(- f_T / 2) ) )
}
# for computing Normalizing constants of G-Wishart distribution according to ATAY-KAYIS AND MASSAM (2005)
I.g = function(G, b = 3, D = diag(ncol(G)), mc = 100)
{
   if (b <= 2) stop("In G-Wishart distribution parameter 'b' has to be more than 2")
   
   G[lower.tri(G, diag = TRUE)] <- 0
   
   sumrowG <- rowSums(G)
   sumcolG <- colSums(G)
   p       <- nrow(G)
   
   Ti      <- chol(solve(D))
   H       <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
   Exp.f_T <- Exp.MC(G, b, H, mc, p)
   
   sumG    <- sum(G)
   c_dT    <- (sumG / 2) * log(pi) + (p * b / 2 + sumG) * log(2) +
              sum(lgamma((b + sumrowG) / 2)) + sum((b + sumrowG + sumcolG) * log(diag(Ti)))
			  
   Ig      <- exp(c_dT) * Exp.f_T
   
   if (is.finite(Ig)){
       return(Ig)
   } else {
       logIg <- c_dT + log (Exp.f_T)
       cat(paste(""), fill = TRUE)
       cat(paste("Normalizing constant is infinte"), fill = TRUE)
	   cat(paste("Log of normalizing constant =", logIg), fill = TRUE)
       cat(paste(""), fill = TRUE)
	   
	   return(logIg) 
   }
}
# sampling from precision matrix K for our BDMCMC algorithm
# according to accept-reject algorithm
sampleK = function(A, b, H, Ts, p, iter = 1)
{
   psi  <- Psi(A, b, H, p)
   cont <- 0
   
   for (i in 1 : (iter * 1000))
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
	
	si <- psi %*% Ts
    return(t(si) %*% si)
}
# Sampling for G-Wishart distribution according to Blocked Gibbs sampling, Wang (2012)
# pairwise block gibbs sampler algorithm for low-dimensional graphs 
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
# pairwise block gibbs sampler algorithm for high-dimensional graphs
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
# computing birth and death rates
ratesf = function(K, A, i, j, b, Ds, Ti, H, method, mc, p)
{
  # (i,j) = 0
  K0       <- K
  K0[i, j] <- 0
  K0[j, i] <- 0     
  K_12     <- K0[j, - j, drop = FALSE]
  K0_ij    <- diag(c(K[i, i], K_12 %*% solve(K0[- j, - j]) %*% t(K_12))) 
  
  # (i,j) = 1
  e     <- c(i, j)
  K_12  <- K[e, - e]  
  K1_ij <- K_12 %*% solve(K[- e, - e]) %*% t(K_12) 

  if (method == "mc")
  {
	Epsi        <- Exp.MC(A = A, b = b, H = H, mc = mc, p = p)
	Ami         <- A
	Ami[i,j]    <- 0
	Epsimi      <- Exp.MC(A = Ami, b = b, H = H, mc = mc, p = p)
	log.ratio.E <- log(Epsi) - log(Epsimi)
  } else {
	   log.ratio.E <- 0
	}
  
  a11    <- K[i, i] - K1_ij[1, 1]
  nustar <- sum(A[i, ])
  
  rate <- sqrt(2) * Ti[i,i] * Ti[j,j] *
		  exp(lgamma((b + nustar) / 2) - lgamma((b + nustar - 1) / 2) 
		  + (1 / 2) * (log(Ds[j, j]) - log(a11))
		  + ((Ds[i, i] - Ds[i, j] ^ 2 / Ds[j, j]) * a11) / 2
		  - sum(Ds[e, e] * (K0_ij - K1_ij)) / 2 + log.ratio.E)
		  
  return(rate)
}
## Main function: BDMCMC algorithm for selecting the best graphs
bdgraph = function(data, n = NULL, npn = "normal", mean = NULL, method = NULL, 
            g.prior = "Uniform", iter = 5000, b = 3, burnin = floor(iter / 2), 
			thin = 1, lambda = NULL, D = NULL, g.start = "full", K.start = NULL, 
			mc = 10, trace = TRUE, save.all = FALSE)
{
  start.time <- Sys.time()
  if (class(data) == "simulate") data <- data $ data
  
  if (is.matrix(data) == FALSE & is.data.frame(data) == FALSE) stop("Data should be a matrix or dataframe")
  if (is.data.frame(data) == TRUE) data <- data.matrix(data)
  if (any(is.na(data))) stop("Data should contain no missing data") 
  if (iter <= burnin)   stop("Number of iteration must be more than number of burn-in")
  
  if (npn != "normal") data <- bdgraph.npn(data = data, npn = npn, npn.thresh = NULL)
  
  dimd <- dim(data)
  
  if (dimd[1] == dimd[2] && sum(abs(data - t(data))) < 1e-5)
  {
     if (is.null(n)) stop("Please specify the number of observations 'n'")
     S <- data
  } else {
     n <- dimd[1]
	 S <- if (!is.null(mean) && mean == 0) t(data) %*% data else n * cov(data)
  }
  
  p <- dimd[2]
  
  gd <- pmatch(g.prior, c("Uniform", "Poisson"))[1]
  
  if(is.na(gd)) gd <- 1
  if(gd == 2 & is.null(lambda)) stop("You should determine value of lambda as the rate of Poisson")

  if (is.null(D))
  { 
     D  <- diag(p)
	 Ti <- D
	 H  <- D
  } else {
	 Ti <- chol(solve(D))
     H  <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
  }
  
  bstar <- b + n
  Ds    <- D + S
  invDs <- solve(Ds)
  Ts    <- chol(invDs)
  Hs    <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))	
	
  if (class(g.start) == "bdgraph") 
  {
     A <- g.start $ last.G
	 K <- g.start $ last.K
  } else {
     A.start <- pmatch(g.start, c("full", "empty", "glasso", "mb", "ct"))[1]	
     if (!is.na(A.start))
	 {
       if (A.start == 1)
	   {
          A               <- 0 * S
          A[upper.tri(A)] <- 1
		  K               <- matrix(rWishart(1, df = bstar + (p - 2), Sigma = invDs), p, p)
	   }
	   
       if (A.start == 2)
	   {
	      A   <- 0 * S
		  psi <- diag(sqrt(rchisq(p, bstar)))
		  si  <- psi %*% Ts
		  K   <- t(si) %*% (si)
		  
		  if (p < 11){
	        K <- block.gibbs.low(K = K, A = A, bstar = bstar, Ds = Ds, p = p)
	      } else { 
	        K <- block.gibbs.high(K = K, A = A, bstar = bstar, Ds = Ds, p = p)
	      }  
	   }
	   
       if (A.start == 3 | A.start == 4 | A.start == 5)
	   {
		  A               <- huge(data, method = g.start)
		  A               <- huge.select(A)
		  A               <- as.matrix(A $ refit)
		  A[lower.tri(A)] <- 0
		  
		  K               <- sampleK(A = A, b = bstar, H = Hs, Ts = Ts, p = p, iter = 1)
	   }
    } else {
	      A                         <- as.matrix(g.start)
          A[lower.tri(A, diag = T)] <- 0
          if (!is.null(K.start)) K <- K.start else K  <- sampleK(A = A, b = bstar, H = Hs, Ts = Ts, p = p, iter = 1) 		  
	}
  }

  if (is.null(method)) method <- ifelse (p < 8, "mc", "fast")

  sample.G   <- all.G       <- vector() # vector of numbers like "10100"
  weights    <- all.weights <- vector() # waiting time for every state
  sumK       <- 0 * K
  
  for (g in 1 : iter)
  {
    if (trace == TRUE && g %% 10 == 0)
	{
	   mes <- paste(c(" iteration: ", g, " from ", iter, ". Graph size= ", sum(A)), collapse="")
   	   cat(mes, "\r")
       flush.console()	
	}
	   
    rates <- 0 * K	
	
    for (i in 1 : (p - 1))
	{
      for (j in (i + 1) : p)
	  {
         if (A[i,j] == 1)
		 {	
		  rates[i,j] <- ratesf(K = K, A = A, i = i, j = j, b = b, Ds = Ds, Ti = Ti, H = H, method = method, mc = mc, p = p)
		  
		  if (gd == 2) rates[i,j] <- (sum(A) / lambda) * rates[i,j]
		 }
		 
		if (A[i,j] == 0)
		{
          Apl       <- A	
          Apl[i, j] <- 1
          Kpl       <- K
		  e         <- c(i, j)
		  
          K12       <- K[e, - e]  
          F         <- K12 %*% solve(K[- e, - e]) %*% t(K12) 

		  a         <- K[e, e] - F

		  sig       <- sqrt(a[1, 1] / Ds[j, j])
		  mu        <- - (Ds[i, j] * a[1, 1]) / Ds[j, j]
          u         <- rnorm(1, mu, sig)
		  v         <- rgamma(1, shape = bstar / 2, scale = 2 / Ds[j, j])

          Kpl       <- K		  
		  Kpl[i, j] <- Kpl[j, i] <- u + F[1, 2]
		  Kpl[j, j] <- v + u ^ 2 / a[1, 1] + F[2, 2]
		  
		  rates[i, j] <- ratesf(K = Kpl, A = Apl, i = i, j = j, b = b, Ds = Ds, Ti = Ti, H = H, method = method, mc = mc, p = p)
		  rates[i, j] <- 1 / rates[i, j]
		  
		  if (gd == 2) rates[i,j] <- (lambda / (sum(A) + 1)) * rates[i,j]
		 }
      }
    }
	
	if (save.all == TRUE)
	{
       indA        <- paste(A[upper.tri(A)], collapse = '')
	   all.G       <- c(all.G, indA)
	   all.weights <- c(all.weights, 1 / sum(rates))
    }
	
    if (g > burnin && g %% thin == 0)
	{
	   sumK <- sumK + K
	   indA <- paste(A[upper.tri(A)], collapse = '')
       wh   <- which(sample.G == indA)
       if (length(weights) != 0 & length(wh) != 0)
	   {
          weights[wh] <- weights[wh] + 1 / sum(rates)
       } else {
          sample.G <- c(sample.G, indA)
          weights  <- c(weights, 1 / sum(rates))
       }
    }
	
    melt      <- cbind(as.matrix(which(upper.tri(rates), arr.ind = TRUE)), rates[upper.tri(rates)])
    rows      <- which(rmultinom(1, 1, melt[,3]) == 1)
    ii        <- melt[rows, 1]
    jj        <- melt[rows, 2]
    A[ii, jj] <- A[ii, jj] + (- 1) ^ (A[ii, jj])
     
	if (p < 11)
	{
	   K <- block.gibbs.low(K = K, A = A, bstar = bstar, Ds = Ds, p = p)
	} else { 
	   K <- block.gibbs.high(K = K, A = A, bstar = bstar, Ds = Ds, p = p)
	}
  }
  
  if (trace == TRUE)
  {
	 mes <- paste(c(" ", iter," iteration done.                    "), collapse = "")
     cat(mes, "\r")
     cat("\n")
     flush.console()
     print(Sys.time() - start.time)  
  }
  
  if (save.all == TRUE)
  {
	 outbdgraph <- list(sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), all.G = all.G, 
		               all.weights = all.weights, last.G = A, last.K = K)
  } else {
	 outbdgraph <- list(sample.G = sample.G, weights = weights, Khat = sumK / (iter - burnin), last.G = A, last.K = K)
  }
  
  class(outbdgraph) <- "bdgraph"
  return(outbdgraph)   
}
# function for computing probability of all links in graph
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
# function for checking the convergency of the BDMCMC algorithm
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
  
  matplot(x = thin * ((1 : length.allG.new)), y = t(ff), type = "l", lty = 1, col = 1,
          xlab = "iteration", ylab = "posterior link probability")
		  
  if (is.null(main)) main <- "Trace plot"
  title(main = main)
  abline(v = thin * length.allG.new / 2, col = "blue")
  
  par(op)
}
# plot size of the graphs for checking the convergency of BDMCMC algorithm
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
	
	if (acf == FALSE & pacf == TRUE){
	   op <- par(mfrow = c(1, 2), pty = "s") 
	   plot(x = 1 : length(all.G), y, type = "l", main = main,
		       ylab = "graph size", xlab = "iteration", ...)
	   abline(h = lin, col = "red")	  
	   pacf(y, main = "PAIC for graph size")
	   par(op)
	}		
}  
# for selecting the most highest posterior probability of the graphs according to bdmcmc result
select = function (output, vis = FALSE)
{
  sample.G <- output $ sample.G
  weights  <- output $ weights
  p        <- nrow(output $ last.G)
  prob.G   <- weights / sum(weights)
  gv       <- c(rep(0, p * (p - 1) / 2))  
  gi       <- sample.G[which(prob.G == max(prob.G))]
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
prob = function(output, g = 2, G = NULL)
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
# for comparing the result with different packages or approaches
roc = function (true, estimate) 
{
	true     <- as.matrix(true)     # true is the adjacency matrix of true graph 
	estimate <- as.matrix(estimate) # estimate is the adjacency matrix of estimated graph 
	true[lower.tri(true, diag = TRUE)]         <- 0
	estimate[lower.tri(estimate, diag = TRUE)] <- 0
	p        <- nrow(true)
	
	tp.all   <- (true != 0) * (estimate != 0) 
	fp.all   <- (true == 0) * (estimate != 0) 
	fn.all   <- (true != 0) * (estimate == 0) 
	tn.all   <- (true == 0) * (estimate == 0)
	
	tp       <- sum(tp.all)
	fp       <- sum(fp.all)
	fn       <- sum(fn.all)
	tn       <- sum(tn.all[upper.tri(tn.all == 1)])
	
	# positive predictive value  
	Precision <- tp / (tp + fp) 
	
	if (is.na(Precision)) Precision <- 0
	# Precision is the probability that a randomly selected link is relevant
	Recall <- tp / (tp + fn) # also called TPR
	
	if (is.na(Recall)) Recall <- 0
	# Recall is the probability that a randomly selected relevant link 
	# is retrieved in a search 
	FPR <- fp / (fp + tn) # False positive rate
	
	if (is.na(FPR)) FPR <- 0
	Accuracy <- (tp + tn) / (tp + tn + fp + fn)
	
	if (is.na(Accuracy)) Accuracy <- 0
	# Specificity <- tn / (tn + fp) # or 1 - false positive rate 
	# Sensitivity <- tp / (tp + fn) # or true positive rate
	# F1score <- 2 * (Precision * Recall) / (Precision + Recall)
	F1score <- (2 * tp) / (2 * tp + fp + fn)
	if (is.na(F1score)) F1score <- 0
	# harmonic mean of precision and recall, called F-measure or balanced F-score:
	roc.matrix <- matrix(c(tp, tn, fp, fn, Recall, FPR, Accuracy, F1score, Precision), 9, 1)
	
	return( round(roc.matrix, 3) )
}
# for comparing the result with different packages or approaches according to the true graph
compare = function (true, estimate, estimate2 = NULL, colnames = NULL, vis = FALSE) 
{
	if (class(true)      == "simulate") true      <- true $ G 
	if (class(true)      == "bdgraph")  {es  <- select(true) ; true <- estimate $ G ; estimate <- es}
	if (class(estimate)  == "bdgraph")  estimate  <- select(estimate) 
	if (class(estimate2) == "bdgraph")  estimate2 <- select(estimate2) 
	if (class(estimate)  == "select")   estimate  <- estimate $ refit
	if (class(estimate2) == "select")   estimate2 <- estimate2 $ refit

	compare.true <- roc(true = true, estimate = true)
	compare.g    <- roc(true = true, estimate = estimate)

	if (!is.null(estimate2))
	{
	   compare.g2  <- roc(true = true, estimate = estimate2)
	   compare.all <- cbind(compare.true, compare.g, compare.g2)
	   if (is.null(colnames)) colnames <- c("True graph", "estimate", "estimate2")
	} else {
		compare.all <- cbind(compare.true, compare.g)
		if (is.null(colnames)) colnames <- c("True graph", "estimate")
	}

   if (vis == TRUE)
   {
	 true     <- graph.adjacency(true,     mode = "undirected", diag = FALSE)
	 estimate <- graph.adjacency(estimate, mode = "undirected", diag = FALSE)
	 
     if (is.null(estimate2))
	 {
	    op <- par(mfrow = c(1, 2), pty = "s")
        plot.igraph(as.matrix(true), layout = layout.circle, main = colnames[1])
		plot.igraph(as.matrix(estimate), layout = layout.circle, main = colnames[2])
	 } else {
	    op        <- par(mfrow = c(2, 2), pty = "s")
		estimate2 <- graph.adjacency(as.matrix(estimate2), mode = "undirected", diag = FALSE)
        plot.igraph(true,      layout = layout.circle, main = colnames[1])
		plot.igraph(estimate,  layout = layout.circle, main = colnames[2])
        plot.igraph(estimate2, layout = layout.circle, main = colnames[3])			
	 }
	 par(op)
   }
   
   colnames(compare.all) <- colnames
   rownames(compare.all) <- c("true positive", "true negative", "false positive", "false negative", 
                "true positive rate", "false positive rate", "accuracy", "balanced F-score", "positive predictive")
				
   return(compare.all)
}
# Data generator according to the graph structure
bdgraph.sim = function(n = 2, p = NULL, graph = NULL, size = NULL, prob = NULL, class = NULL, 
                 v = NULL, u = NULL, G = NULL, K = NULL, sigma = NULL, mean = 0, vis = FALSE)
{
    if (is.null(graph) & is.null(G)) graph <- "random"
	
	if (!is.null(G))
	{
	    graph <- "fixed"
		p     <- nrow(G)
    } 

	if (is.null(p)) p <- 10
	
	if (graph == "random" & is.null(G))
	{
	   G <- matrix(0, p, p)
	   
	   if (is.null(size))
	   {
		  if (is.null(prob)) prob <- 0.2
		  if (prob < 0 | prob > 1) stop("'prob' should be between zero and one")
		  G[upper.tri(G)] <- rbinom(p * (p - 1) / 2, 1, prob)
	    } else {
		  if (size < 0 | size > p * (p - 1) / 2)  stop("Graph size should be between zero and p*(p-1)/2")
          smp <- sample(1 : (p * (p - 1) / 2), size, replace = FALSE)
          G[upper.tri(G)][smp] <- 1
		}
	   G <- G + t(G)
	}
	
	if (graph == "cluster" & is.null(G))
	{
	  # partition variables
	  if (is.null(class))
	  { 
	     if (!is.null(size)){
			class <- length(size)
		 } else {
			class <- max(2, ceiling(p / 20))
		 }
	  }
	  
	  g.large <- p %% class
	  g.small <- class - g.large
	  n.small <- floor(p / class)
	  n.large <- n.small + 1
	  vp      <- c(rep(n.small, g.small), rep(n.large, g.large)) 
	  G       <- matrix(0, p, p)
	  	 
	 if (is.null(size))
	 {
		if (is.null(prob)) prob <- 0.2
		if (prob < 0 | prob > 1) stop("'prob' should be between zero and one")

		for (i in 1 : class)
		{
		   tmp <- if (i == 1) (1 : vp[1]) else ((sum(vp[1 : (i-1)]) + 1) : sum(vp[1:i]))
		   gg                <- matrix(0, vp[i], vp[i])
           gg[upper.tri(gg)] <- rbinom(vp[i] * (vp[i] - 1) / 2, 1, prob)
		   G[tmp, tmp]       <- gg
		}
	 } else {
		
		  if (class != length(size))  stop("Number of graph sizes is not match with number of clusters")
		  if (sum(size) < 0 | sum(size) > p * (p - 1) / 2)   stop("Total graph sizes should be between zero and p*(p-1)/2")
		  
		  for (i in 1 : class)
		  {
		     tmp <- if (i == 1) (1 : vp[1]) else ((sum(vp[1 : (i-1)]) + 1) : sum(vp[1:i]))
		     gg  <- matrix(0, vp[i], vp[i])
		     smp <- sample(1 : (vp[i] * (vp[i] - 1) / 2), size[i], replace = FALSE)
             gg[upper.tri(gg)][smp] <- 1
			 G[tmp, tmp]            <- gg
		  }
		}
	   G <- G + t(G)	   
	}
	
	if (graph == "circle" & is.null(G))
	{
	    G       <- toeplitz(c(0, 1, rep(0, p - 2)))
        G[1, p] <- 1
		G[p, 1] <- 1
	}
	
    if (is.null(sigma) & is.null(K))
	{
        if(is.null(u)) u <- 0.2
	    if(is.null(v)) v <- 0.3	  	
		K       <- G * v
		# make K positive definite and standardized
		diag(K) <- abs(min(eigen(K) $ values)) + u
		sigma   <- cov2cor(solve(K))
		K       <- solve(sigma)
	} else {
	    if (is.null(sigma)) sigma <- solve(K)
		if (is.null(K))     K     <- solve(sigma)
		G <- 1 * (abs(K) > 0.02)
		p <- nrow(G)
	}
	diag(G) <- 0
	
	# generate multivariate normal data
	if (typeof(mean) == "double") mean <- rep(mean, p)
	R <- chol(sigma)
	z <- matrix(rnorm(p * n), p, n)
	d <- t(R) %*% z + mean
	
	# graph visualization
	if (vis == TRUE)
	{
	  graphG <- graph.adjacency(G, mode = "undirected", diag = FALSE)
	  if (p < 20){
        plot.igraph(graphG, layout = layout.circle, main = "True graph structure",
		            edge.color = 'black', vertex.color = "white", vertex.size = 10)
	  } else {
	    plot.igraph(graphG, layout = layout.circle, main = "True graph structure", 
	                edge.color = 'black', vertex.color = "white", vertex.size = 2)
	  }
	}
	
	dimlab     <- as.character(1 : p)
	simulation <- list(data = t(d), sigma = sigma, K = K, G = Matrix(G, sparse = TRUE, dimnames = list(dimlab, dimlab)), graph = graph)
	
	class(simulation) <- "simulate"
	return(simulation)
}
# print function of simulation data
print.simulate = function (x, ...)
{
  p <- ncol(x $ sigma)
  
  cat(paste("  Data generated by bdgraph.sim"), fill = TRUE)
  cat(paste("  Sample size     =", nrow(x $ data)), fill = TRUE)
  cat(paste("  Graph type      =", x $ graph), fill = TRUE)
  cat(paste("  Number of nodes =", p), fill = TRUE)
  cat(paste("  Graph size      =", sum(x $ G) / 2), fill = TRUE)
  cat(paste("  Sparcity        =", round(sum(x $ G) / (p * (p - 1)), 4)), fill = TRUE)
}
# plot for class "simulate" from bdgraph.sim function
plot.simulate = function(x, main = NULL, layout = layout.circle, ...)
{
    if (is.null(main)) main <- "True graph structure"
  	g <- graph.adjacency(as.matrix(x $ G), mode = "undirected", diag = FALSE)
	
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
    plot.igraph(G, layout = layout, main = main, 
     sub = paste(c("Posterior probability = ", round(sort(prob.G, decreasing = TRUE)[i], 4)), collapse = ""), ...)	   
  }
  if (g > 1 & g < 7) par(op)
}
# summary of the result according to bdgraph
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
# print of the result according to bdgraph
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
rGWishart = function(n = 1, b = 3, D = diag(2), G = NULL, method = "block gibbs", start.K = NULL)
{
  if (sum((G == 1) * (G == 0)) != 0) stop("Elements of matrix G should be zero or one")
  if (!is.null(G)) G[lower.tri(G, diag(TRUE))] <- 0
  p <- nrow(D)
  
  if (is.null(G)) 
  {
      method <- "rWishart"
  } else {
      if (sum(upper.tri(G)) == sum(G[upper.tri(G == 1)])) method <- "rWishart"
  }	
  
  id      <- pmatch(method, c("block gibbs", "accept-reject", "rWishart"))[1]
  samples <- array(0, c(p, p, n))
  
  if (id == 1) 
  {
    if (is.null(start.K))
	{
	   Ts      <- chol(solve(D))
       H       <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))
	   start.K <- sampleK(A = G, b = b, H = H, Ts = Ts, p = p, iter = 1)
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
  
  if (id == 2) 
  {
    Ts <- chol(solve(D))
    H  <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))
	
    for (i in 1 : n)
	{
       samples[,,i] <- sampleK(A = G, b = b, H = H, Ts = Ts, p = p, iter = 1)
    }
  }
  
  if (id == 3) samples <- rWishart(n = n, df = b + p - 1, Sigma = solve(D))
  
  return(samples)   
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
		data <- qnorm(pmin(pmax(apply(data, 2, rank) / n, npn.thresh), 1 - npn.thresh))
    	data <- data / sd(data[ , 1])
	}

	if(npn == "skeptic") data <- 2 * sin(pi / 6 * cor(data, method = "spearman"))
	
	return(data)
}






