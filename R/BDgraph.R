## All functions for the "BDgraph" package
###########################################################
## function for only obtaining elements of Psi matrix
Psi = function(A, b, H, p)
{
  nu <- apply(A, 1, sum)
  psi <- diag(sqrt(rchisq(p, b + nu)))
  psi[A == 1] <- rnorm(1)
  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
        if (A[i, j] == 0){
            psi[i, j] <- - sum(psi[i, i : (j - 1)] * H[i : (j - 1), j])
            if (i > 1){
               for (r in 1 : (i - 1)){
                 psi[i, j] <- psi[i, j] - ((sum(psi[r, r : i] * H[r : i, i])) *
                             (sum(psi[r, r : j] * H[r : j, j]))) / (psi[i, i])
               }
            }
        }
    }
  }
  return(psi)
}
# Algorithm 3.1: for Monte Carlo approxiation for expection in normalizing constant
Exp.MC = function(A, b, H, MCiter, p)
{
   f_T <- vector()
   for (i in 1 : MCiter){
       psi <- Psi(A = A, b = p, H = H, p = p)
	   dd <- A
	   dd[lower.tri(dd == 0, diag = T)] <- 1
       f_T[i] <- exp(- sum((1 - dd) * psi * psi) / 2)
   }
   return(mean(f_T))
}
# for computing Normalizing constans of G-Wishart distribution according to ATAY-KAYIS AND mASSAM (2005)
I.g = function(A, b, D, MCiter = 500)
{
   if (b <= 2){
      stop("parameter 'b' in G-Wishart distribution has to be more than 2")
   }
   p <- nrow(A)
   Ti <- chol(solve(D))
   H <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
   Exp.f_T <- Exp.MC(A, b, H, MCiter, p)
   c_dT <- 0
   for (i in 1 : p){
      c_dT <- c_dT + ((sum(A[i, ]) / 2) * log(pi) + ((b + 2 * sum(A[i, ])) / 2) * log(2) +
      lgamma((b + sum(A[i, ])) / 2) + (b + sum(A[i, ]) + sum(A[ , i])) * log(Ti[i, i]))
   }
   c_dT <- exp(c_dT)
   cat(paste(""), fill = TRUE)
   cat(paste(c("Normalizing constant = ", c_dT * Exp.f_T),collapse = ""), fill = TRUE)
}
# sampling from precistion matrix K for our BDMCMC algorithm
# according to accept-reject algorithm
sampleK = function(A, b, H, Ts, p, iter = 1)
{
   psi <- Psi(A, b, H, p)
   cont <- 0
   for (i in 1 : (iter * 1000)){  
      psi.new <- Psi(A, b, H, p)
	  alpha <- exp((sum((1 - (diag(p) + A)) * psi * psi) - 
	  sum((1 - (diag(p) + A)) * psi.new * psi.new)) / 2)
	  if (is.nan(alpha)) alpha <- 0
	  if (runif(1) < alpha){ 
	    cont <- cont + 1
		psi.new <- psi
	    }  
      if (cont == iter) break
	}
    return(t(psi %*% Ts) %*% (psi %*% Ts))
}
# Sampling for G-Wishart distribution according to Blocked Gibbs sampling, Wang (2012)
block.gibbs = function (K, A, bstar, Ds, p, gibbs.iter = 1)
{
for (k in 1 : gibbs.iter){
  for (i in 1 : (p - 1)){
	if (sum(A[i, ]) != 0){
      for (j in (i + 1) : p){
	    if (A[i, j] == 1){
		  B <- Ds[c(i, j), c(i, j)]
		  a <- rWishart(1, df = bstar + (p - 2), Sigma = solve(B))
		  Kvv <- K[- c(i, j), - c(i, j)]
          kjv <- K[c(i, j), - c(i, j)]
          Kc <- matrix(a, 2, 2) + (kjv) %*% (solve(Kvv)) %*% t(kjv)
          K[c(i, j), c(i, j)] <- (Kc + t(Kc)) / 2 # Numerical stable
		}
      }
	}
    if (sum(A[i, ]) + sum(A[ , i]) == 0){
	  a <- rgamma(1, (bstar + (p - 1)) / 2, Ds[i, i] / 2)
	  Kvv <- K[- i, - i]
	  kjv <- K[i, - i, drop = FALSE]
      K[i, i] <- a + (kjv) %*% (solve(Kvv)) %*% t(kjv)
	}
  }
  if (sum(A[ , p]) == 0){
	a <- rgamma(1, (bstar + (p - 1)) / 2, Ds[p, p] / 2)
	Kvv <- K[- p, - p]
    kjv <- K[p, - p, drop = FALSE]
    K[p, p] <- a + (kjv) %*% (solve(Kvv)) %*% t(kjv)
  }
}
return(K)
}
#auxiliary function to get covariance matrix
get.S = function(data, n, tol = 1e-5, meanzero)
{
  if (ncol(data) != nrow(data)){
     n <- nrow(data)
	 if (meanzero == TRUE) S <- t(data) %*% data
     if (meanzero == FALSE) S <- n * cov(data)
  } else {
     if (sum(abs(data - t(data))) > tol){
        n <- nrow(data)
	    if (meanzero == TRUE) S <- t(data) %*% data
        if (meanzero == FALSE) S <- n * cov(data)
    }
  }
  return(list(S = S, n = n))
}
# Algorithm 2.1: BD-MCMC algorithm for low-dimentional problem (roughly graphs with more less 8 nodes)
bdgraph.low = function(data, n = NULL, meanzero = FALSE, iter = 5000, burnin = floor(iter / 2), skip = 1, gamma.b = 1, 
prior.g = "Uniform", b = 3, D = NULL, start.g = "full", MCiter = 10, summary = FALSE, verbose = TRUE, save.all = FALSE, last.objects = NULL, time = TRUE)
{
  start.time <- Sys.time()
  if (is.matrix(data) == F & is.data.frame(data) == F){
     stop("Data should be a matrix or dataframe")
  }
  if (is.data.frame(data) == T) data <- data.matrix(data)
  if (any(is.na(data))) stop("Data should contain no missing data") 
  if (iter <= burnin){
    stop("Number of iterations have to be more than the number of burn-in iterations")
  }
  if (gamma.b <= 0){
    stop("birth rate 'gamma.b' has to be positive value")
  }
  Sn <- get.S(data = data, n = n, meanzero = meanzero)
  if (is.null(Sn $ n) & is.null(n)){
    stop("You have to specify the number of observations 'n'")
  }
  S <- Sn $ S
  n <- Sn $ n
  p <- ncol(S)
  id <- pmatch(prior.g, c("Uniform", "Poisson"))[1]
  if(!is.na(id)){
     if(id == 1) pr <- 0
     if(id == 2) pr <- 1
    }
  if (is.null(last.objects)){
    start.A <- pmatch(start.g, c("full", "empty", "glasso"))[1]	
    if (!is.na(start.A)){
       if (start.A == 1){
          A <- 0 * S
          A[upper.tri(A)] <- 1
	      }
       if (start.A == 2) A <- 0 * S
       if (start.A == 3){
		  A <- huge(data)
		  A <- huge.select(A)
		  A <- as.matrix(Matrix(A $ refit, sparse = F))
		  A[lower.tri(A)] <- 0
		}
    }
  } else {
       A <- last.objects $ A
  }
  if (is.null(D)) D <- diag(p)
  bstar <- b + n
  Ds <- D + S
  invDs <- solve(Ds)
  Ts <- chol(invDs)
  Ti <- chol(D)
  H <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
  Hs <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))
  if (is.null(last.objects)){
    if (sum(A[upper.tri(A)]) == p * (p - 1) / 2){
       K <- matrix(rWishart(1, df = bstar + (p - 2), Sigma = invDs), p, p)
       } else {
       K <- sampleK(A, bstar, Hs, Ts, p, iter = 1)
       }
    } else {
       K <- last.objects $ K
    }
  As <- allA  <- vector() # vector of numbers like "10100"
  lambda <- all.lambda <- vector() # waiting time for every state
  alla <- ceiling(iter / 2000)# for saving save.all which we need it for plotcoda function
  sumK <- 0 * K
  for (g in 1:iter){
    if(verbose == T){
	   mes <- paste(c("    MCMC iterations : ", g, " from ", iter, ". Graph size = ", sum(A)), collapse="")
   	   cat(mes, "\r")
       flush.console()	
	}      
    rates <- 0 * K
    for (i in 1:(p-1)){
     for (j in (i+1):p){
       if (A[i,j] == 0) rates[i,j] <- gamma.b
       if (A[i,j] == 1){
        mu <- - (Ds[i,j] * K[i,i]) / Ds[j,j]
	    sig <- sqrt(K[i,i] / Ds[j,j])
        k_xi <- rnorm(1, mu, sig)
        b_xi <- dnorm(k_xi, mu, sig)
        nustar <- sum(A[i,])
        Epsi <- Exp.MC(A = A, b = b, H = H, MCiter = MCiter, p = p)
        Ami <- A
        Ami[i,j] <- 0
        Epsimi <- Exp.MC(A = Ami, b = b, H = H, MCiter = MCiter, p = p)
		Kmi <- K
		Kmi[i,j] <- Kmi[j,i] <- 0
		Kmi[j,j] <- K[j,j] + get.cs_c(K = K, p = p, i = i, j = j)
        if (sum(A) == 0 & pr == 0) pr <- 1
        rates[i, j] <- (((sum(A)) ^ pr) * ((gamma.b) ^ (1 - pr))) * 
		(b_xi) * 2 * sqrt(pi) * Ti[i,i] * Ti[j,j] *
        exp(lgamma((bstar + nustar) / 2) - lgamma((bstar + nustar - 1) / 2) +
        log(Epsi) - log(Epsimi) + ((bstar - 2) / 2) * (log(det(Kmi)) - log(det(K))) -
        (sum(Ds * t(Kmi - K))) / 2) 
        if (is.finite(rates[i,j]) == FALSE) rates[i,j] <- gamma(170)
        }
     }
    }
	if ( save.all == TRUE & g %% alla == 0){
        indA <- paste(A[upper.tri(A)], collapse = '')
		allA <- c(allA, indA)
		all.lambda <- c(all.lambda, 1 / sum(rates))
    }
    if (g > burnin && g %% skip == 0){
	  sumK <- sumK + K
	  indA <- paste(A[upper.tri(A)], collapse = '')
      wh <- which(As == indA)
      if (length(lambda) != 0 & length(wh) != 0){
         lambda[wh] <- lambda[wh] + 1 / sum(rates)
      } else {
         As <- c(As, indA)
         lambda <- c(lambda, 1 / sum(rates))
        }
    }
	melt <- cbind(as.matrix(which(upper.tri(rates), arr.ind = TRUE)), rates[upper.tri(rates)])
    rows <- which(rmultinom(1, 1, melt[ , 3]) == 1)
    ii <- melt[rows, 1]
    jj <- melt[rows, 2]
    A[ii,jj] <- A[ii,jj] + (- 1) ^ (A[ii,jj])  
    K <- block.gibbs(K = K, A = A, bstar = bstar, Ds = Ds, p, gibbs.iter = 1)
  }
  if(verbose == TRUE){
	mes = paste(c("    ", iter," iterations done.                              "), collapse = "")
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  if(time == TRUE) print(Sys.time() - start.time)  
  last.objects <- list(K = K, A = A)  
  if (summary == FALSE){
     if (save.all == TRUE){
	    bdgraphres <- c(list(As = As, lambda = lambda, Khat = sumK / (iter - burnin), p = p, allA = allA, 
		              all.lambda = all.lambda, alla = alla), list(last.objects = last.objects))
	    class(bdgraphres) <- "bdgraph"
        return(bdgraphres)
	 } else {
	    bdgraphres <- c(list(As = As, lambda = lambda, Khat = sumK / (iter - burnin), p = p), 
		                list(last.objects = last.objects))
	    class(bdgraphres) <- "bdgraph"
        return(bdgraphres)
	 }
  } else {
        output <- list(As = As, lambda = lambda, p = p)
		class(output) <- "bdgraph"
		print.bdgraph(output)
        # for phat
        pvec <- c(rep(0, p * (p - 1) / 2))
        for (i in 1:length(As)){
           inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
	       pvec[inp] <- pvec[inp] + lambda[i]
        }
        phat <- matrix(0, p, p)
        phat[upper.tri(phat)] <- pvec / sum(lambda)
		Khat <- sumK / ((iter - burnin) * skip)
		return.list <- list(round(phat, 2), round(Khat, 3))
		names(return.list) <- c("Posterior edge inclusion ","Estimation for precision matrix")
        return(return.list)	
    }
}
# function for computing "c.star - c" from the matrix K for k_jj
get.cs_c = function(K, p, i, j){
   kjj <- K[c((1:p)[- j]), c((1:p)[- j])]
   kjv <- K[j,][- j]
   B <- solve(kjj)
   return(K[i,j] * (K[i,j] * B[i,i] - 2 * (kjv %*% B[,i])))  
}
## Algorithm 2.1: BD-MCMC algorithm for high-dimentional problem (roughly graphs with more than 8 nodes)
bdgraph.high = function(data, n = NULL, meanzero = FALSE, iter = 5000, burnin = floor(iter / 2), skip = 1, 
              gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, start.g = "full", summary = FALSE, verbose = TRUE, 
			  save.all = FALSE, last.objects = NULL, time = TRUE)
{
  start.time <- Sys.time()
  if (is.matrix(data) == F & is.data.frame(data) == F){
     stop("Data should be a matrix or dataframe")
  }
  if (is.data.frame(data) == T) data <- data.matrix(data)
  if (any(is.na(data))) stop("Data should contain no missing data") 
  if (iter <= burnin){
    stop("Number of iterations have to be more than the number of burn-in iterations")
  }
  if (gamma.b <= 0){
    stop("birth rate 'gamma.b' has to be positive value")
  }
  Sn <- get.S(data = data, n = n, meanzero = meanzero)
  if (is.null(Sn $ n) & is.null(n)){
    stop("If you provide the covariance matrix, you should specify the number of observations")
  }
  S <- Sn $ S
  n <- Sn $ n
  p <- ncol(S)
  id <- pmatch(prior.g, c("Uniform", "Poisson"))[1]
  if (!is.na(id)){
     if (id == 1) pr <- 0
     if (id == 2) pr <- 1
    }
  if (is.null(last.objects)){
    start.A <- pmatch(start.g, c("full", "empty", "glasso"))[1]	
    if (!is.na(start.A)){
       if (start.A == 1){
          A <- 0 * S
          A[upper.tri(A)] <- 1
	      }
       if (start.A == 2) A <- 0 * S
       if (start.A == 3){
		  A <- huge(data)
		  A <- huge.select(A)
		  A <- as.matrix(Matrix(A $ refit, sparse = F))
		  A[lower.tri(A)] <- 0
		}
    }
  } else {
       A <- last.objects $ A
  }
  if (is.null(D)) D <- diag(p)
  Ti <- chol(solve(D))
  bstar <- b + n
  Ds <- D + S
  invDs <- solve(Ds)
  Ts <- chol(invDs)
  Hs <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))
  # for starting point of precision matrix K
  if (is.null(last.objects)){
    if (sum(A[upper.tri(A)]) == p * (p - 1) / 2){
       K <- matrix(rWishart(1, df = bstar + (p - 2), Sigma = invDs), p, p)
       } else {
       K <- sampleK(A, bstar, Hs, Ts, p, iter = 1)
       }
    } else {
       K <- last.objects $ K
    }
  As <- allA  <- vector() # vector of numbers like "10100"
  lambda <- all.lambda <- vector() # waiting time for every state
  alla <- ceiling(iter / 2000) # for saving allA which we need it for plotcoda function
  sumK <- 0 * K
  for (g in 1:iter){
    if(verbose == T){
	   mes <- paste(c("    MCMC iterations : ", g, " from ", iter, ". Graph size = ", sum(A)), collapse="")
   	   cat(mes, "\r")
       flush.console()	
	}     
    rates <- 0 * K
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        if (A[i,j] == 0) rates[i,j] <- gamma.b
        if (A[i,j] == 1){
		  sig <- sqrt(K[i,i] / Ds[j,j])
          mu <- - (Ds[i,j] * K[i,i]) / Ds[j,j]
          k_xi <- rnorm(1, mu, sig)
          b_xi <- dnorm(k_xi, mu, sig)
	      Kmi <- K
          Kmi[i,j] <- Kmi[j,i] <- 0
		  Kmi[j,j] <- K[j,j] + get.cs_c(K = K, p = p, i = i, j = j) 
          nustar <- sum(A[i,])
		  if (sum(A) == 0 & pr == 0) pr <- 1	 
		  rates[i,j] <- (((sum(A)) ^ pr) * ((gamma.b) ^ (1 - pr))) * 
		  (b_xi) * 2 * sqrt(pi) * Ti[i,i] * Ti[j,j] *
          exp(lgamma((bstar + nustar) / 2) - lgamma((bstar + nustar - 1) / 2) +
		  ((bstar - 2) / 2) * (log(det(Kmi)) - log(det(K))) -
                     (sum(Ds * t(Kmi - K))) / 2 )
          if (is.finite(rates[i,j]) == FALSE) rates[i,j] <- gamma(170)
        }
      }
    }
	if ( save.all == TRUE & g %% alla == 0){
        indA <- paste(A[upper.tri(A)], collapse = '')
		allA <- c(allA, indA)
		all.lambda <- c(all.lambda, 1 / sum(rates))
    }
    if (g > burnin && g %% skip == 0){
	  sumK <- sumK + K
	  indA <- paste(A[upper.tri(A)], collapse = '')
      wh <- which(As == indA)
      if (length(lambda) != 0 & length(wh) != 0){
         lambda[wh] <- lambda[wh] + 1 / sum(rates)
      } else {
         As <- c(As, indA)
         lambda <- c(lambda, 1 / sum(rates))
        }
    }
    melt <- cbind(as.matrix(which(upper.tri(rates), arr.ind = TRUE)), rates[upper.tri(rates)])
    rows <- which(rmultinom(1, 1, melt[,3]) == 1)
    ii <- melt[rows,1]
    jj <- melt[rows,2]
    A[ii,jj] <- A[ii,jj] + (- 1) ^ (A[ii,jj]) 
	K <- block.gibbs(K, A, bstar, Ds, p, gibbs.iter = 1)
  }
  if(verbose == TRUE){
	mes = paste(c("    ", iter," iterations done.                              "), collapse = "")
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }
  if(time == TRUE) print(Sys.time() - start.time)
  last.objects <- list(K = K, A = A)  
  if (summary == FALSE){
     if (save.all == TRUE){
	    bdgraphres <- c(list(As = As, lambda = lambda, Khat = sumK / (iter - burnin), p = p, allA = allA, 
		              all.lambda = all.lambda, alla = alla), list(last.objects = last.objects))
	    class(bdgraphres) <- "bdgraph"
        return(bdgraphres)
	 } else {
	    bdgraphres <- c(list(As = As, lambda = lambda, Khat = sumK / (iter - burnin), p = p), 
		                list(last.objects = last.objects))
	    class(bdgraphres) <- "bdgraph"
        return(bdgraphres)
	 }
  } else {
        output <- list(As = As, lambda = lambda, p = p)
		class(output) <- "bdgraph"
		print.bdgraph(output)
        # for phat
        pvec <- c(rep(0, p * (p - 1) / 2))
        for (i in 1:length(As)){
           inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
	       pvec[inp] <- pvec[inp] + lambda[i]
        }
        phat <- matrix(0, p, p)
        phat[upper.tri(phat)] <- pvec / sum(lambda)
		Khat <- sumK / ((iter - burnin) * skip)
		return.list <- list(round(phat, 2), round(Khat, 3))
		names(return.list) <- c("Posterior edge inclusion ", "Estimation for precision matrix")
        return(return.list)	
    }
}
## Main function: BDMCMC algorithm for selecting the best graphs
bdgraph = function(data, n = NULL, meanzero = FALSE, model = NULL, iter = 5000, burnin = floor(iter / 2), skip = 1, gamma.b = 1, 
prior.g = "Uniform", b = 3, D = NULL, start.g = "full", MCiter = 10, summary = FALSE, verbose = TRUE, save.all = FALSE, last.objects = NULL, time = TRUE)
{
  if (is.null(model)){
     p <- ncol(data)
	 if (p < 8) {
		model <- "low"
	 } else {
	    model <- "high"
		}
  }
  if (model == "low"){
    return(bdgraph.low(data, n, meanzero, iter, burnin, skip, gamma.b, prior.g, b, D, start.g, MCiter, summary, verbose, save.all, last.objects, time))
  } else {
    return(bdgraph.high(data, n, meanzero, iter, burnin, skip, gamma.b, prior.g, b, D, start.g, summary, verbose, save.all, last.objects, time))
  }
}
# function for comuting probability of all links in graph
phat = function(output, round = 3)
{
   As <- output $ As
   lambda <- output $ lambda
   p <- output $ p
   pvec <- c(rep(0, p * (p - 1) / 2))
   for (i in 1:length(As)){
      inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
	  pvec[inp] <- pvec[inp] + lambda[i]
   }
   phat <- matrix(0, p, p)
   phat[upper.tri(phat)] <- pvec / sum(lambda)
   cat(paste(""), fill = TRUE)
   cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
   cat(paste(""), fill = TRUE)
   return(round(phat, round))
}
# plot for probability of graphs according to number of their links
densplot = function(output, xlim = NULL, ylim = NULL, main = NULL)
{
   p <- output $ p
   As <- output $ As
   lambda <- output $ lambda
   suma <- sapply(As, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
   x <- unique(suma)
   lambdag <- vector()
   for (i in 1:length(x)){
      lambdag[i] <- sum(lambda[which(suma == x[i])])
   }
   plot(x = x, y = lambdag, type = "h", main = main, xlim = xlim, ylim = ylim,
   ylab = "Pr(number of links in the graph|data)", xlab = "number of links in the graph")
}
# function for checking the convergency of the BDMCMC algorithm
plotcoda = function(output, skip = ceiling(length(output $ allA) / 2000), verbose = TRUE, main = NULL)
{
  p <- output $ p
  all.lambda <- output $ all.lambda
  allA <- output $ allA
  if (skip == 1) skip = 2
  allA.new <- allA[c(skip * (1 : floor(length(allA) / skip)))]
  all.lambda.new <- all.lambda[c(skip * (1 : floor(length(all.lambda) / skip)))]
  length.allA.new <- length(allA.new)
  ff <- matrix(0, p * (p - 1) / 2, length.allA.new)
  inp <- which(unlist(strsplit(as.character(allA.new[1]), "")) == 1)
  ff[inp,1] <- 1
  ffv <- 0 * ff[ , 1]
  ffv[inp] <- all.lambda.new[1]
  for (g in 2 : length.allA.new){
     if (verbose == TRUE){
	   mes <- paste(c("Calculating cumulative occupancy fractions....in progress : ", floor(100 * g / length.allA.new), "%"), collapse = "")
	   cat(mes, "\r")
	   flush.console()	
     }
	 inp <- which(unlist(strsplit(as.character(allA.new[g]), "")) == 1)
	 ffv[inp] <- ffv[inp] + all.lambda.new[g]
	 ff[,g] <- ffv / sum(all.lambda.new[c(1:g)])    	 
    }  
  if(verbose == TRUE){
	mes = paste(c("Calculating cumulative occupancy fractions....done.                   "), collapse = "")
    cat(mes, "\r")
    cat("\n")
    flush.console()
  } 
  matplot(x = skip * ((output $ alla) * (1:length.allA.new)), y = t(ff), type = "l", lty = 1, col = 1,
  xlab = "number of iterations", ylab = "cumulative occupancy fractions for each links")
  if (is.null(main)) main = "To check the convergency of the BDMCMC algorithm"
  title(main = main)
  abline(v = skip * (output $ alla) * length.allA.new / 2, col = "blue")
}
# plot size of the graphs for checking the convergency of BDMCMC algorithm
traceplot = function(output, acf = FALSE, pacf = FALSE, xlim = NULL, ylim = NULL, main = NULL)
{
    allA <- output $ allA
    y <- sapply(allA, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
	if (acf == FALSE & pacf == FALSE){
       plot(x = (output $ alla) * (1:length(allA)), y, type = "l", main = main,
            ylab = "sum of links in the graphs", xlab = "iterations", xlim = xlim, ylim = ylim)	
	}
	if (acf == TRUE & pacf == TRUE){
	   op <- par(mfrow = c(2, 2), pty = "s") 
	   plot(x = (output $ alla) * (1 : length(allA)), y, type = "l", main = main,
		       ylab = "sum of links in the graphs", xlab = "iterations", xlim = xlim, ylim = ylim)
	   acf(y, main = "AIC for graph sizes")
	   pacf(y, main = "PAIC for graph sizes")
	   par(op)
	}		
	if (acf == TRUE & pacf == FALSE){
	   op <- par(mfrow = c(1, 2), pty = "s") 
	   plot(x = (output $ alla) * (1 : length(allA)), y, type = "l", main = main,
		       ylab = "sum of links in the graphs", xlab = "iterations", xlim = xlim, ylim = ylim)
	   acf(y, main = "AIC for graph sizes")
	   par(op)
	}
	if (acf == FALSE & pacf == TRUE){
	   op <- par(mfrow = c(1, 2), pty = "s") 
	   plot(x = (output $ alla) * (1 : length(allA)), y, type = "l", main = main,
		       ylab = "sum of links in the graphs", xlab = "iterations", xlim = xlim, ylim = ylim)
	   pacf(y, main = "PAIC for graph sizes")
	   par(op)
	}		
}  
# for selecting the most highest posterior probability of the graphs according to bdmcmc result
select = function (output, g = 1, mode = "circle", edge.col = "black", label.col = "black", vertex.col = "red")
{
  As <- output $ As
  lambda <- output $ lambda
  p <- output $ p
  prob.A <- lambda / sum(lambda)
  graphi <- list()
  gv <- c(rep(0, p * (p - 1) / 2))  
  for (i in 1 : g){
    dev.new()
    gi <- As[[which(prob.A == sort(prob.A, decreasing = T)[i])]]
    gv <- 0 * gv
    gv[which(unlist(strsplit(as.character(gi), "")) == 1)] = 1
    graphi[[i]] <- matrix(0, p, p)
    graphi[[i]][upper.tri(graphi[[i]])] <- gv
    G <- network(graphi[[i]], directed = F)
    if (i == 1){
       main = paste(c("BEST GRAPH with size ", sum(graphi[[1]]), " : graph with highest probability"), collapse = "") 
    } else {
	   main = paste(c("Graph with ", i, "th highest probability : graph size=", sum(graphi[[i]])), collapse = "")
	   }
    plot.network(G, label = network.vertex.names(G), mode = mode, main = main,
	 edge.col = edge.col, label.col = label.col, vertex.col = vertex.col,
     sub = paste(c("Posterior probability of graph = ", round(sort(prob.A, decreasing = TRUE)[i], 4)), collapse = ""))
  }
  return(Matrix(graphi[[1]] + t(graphi[[1]]), sparse = TRUE))
}
# computing the probability of all the possible graphs or one specific graph 
prob = function(output, g = 2, A = NULL)
{
   As <- output $ As
   lambda <- output $ lambda
   if (is.null(A)){
      p <- output $ p
      prob.A <- lambda / sum(lambda)
      graphi <- list()
      gv <- c(rep(0, p * (p - 1) / 2))        
      for (i in 1 : g){
	    gi <- As[[which(prob.A == sort(prob.A, decreasing = T)[i])]]
        gv <- 0 * gv
        gv[which(unlist(strsplit(as.character(gi), "")) == 1)] = 1
        graphi[[i]] <- matrix(0, p, p)
        graphi[[i]][upper.tri(graphi[[i]])] <- gv
	  }
      return(list(list.A = graphi, prob.A = prob.A[1 : g]))
   } else {
      indA <- paste(A[upper.tri(A)], collapse = '')
      lambda.g <- 0
      for (i in 1:length(As)){
          if (identical(indA,As[i]) == TRUE) lambda.g <- lambda[i]
	    }
      cat(paste(""), fill = TRUE)
      mes = paste(c("  Posterior probability of the graph = ", lambda.g / sum(lambda)), collapse = "")
      cat(mes, "\r")
      cat("\n")
    }
}
# for comparing the result with different packages or approaches
roc = function (true.g, est.g) 
{
	true.g <- as.matrix(true.g) # true.g is the adjacency matrix of true graph 
	est.g <- as.matrix(est.g)   # est.g is the adjacency matrix of estimated graph 
	true.g[lower.tri(true.g, diag = T)] <- 0
	est.g[lower.tri(est.g, diag = T)] <- 0
	p <- nrow(true.g)
	tp.all <- (true.g != 0) * (est.g != 0) 
	fp.all <- (true.g == 0) * (est.g != 0) 
	fn.all <- (true.g != 0) * (est.g == 0) 
	tn.all <- (true.g == 0) * (est.g == 0)
	tp <- sum(tp.all)
	fp <- sum(fp.all)
	fn <- sum(fn.all)
	tn <- sum(tn.all) - (p * (p + 1)) / 2
	# positive predictive value  
	Precision <- tp / (tp + fp) 
	if (is.na(Precision)) Precision <- 0
	# Precision is the probabilty that a randomly selected link is relevant
	Recall <-  tp / (tp + fn) # also called TPR
	if (is.na(Recall)) Recall <- 0
	# Recall is the probability that a randomely selected relevant link 
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
	roc.list <- list(true.positive = tp, true.negative = tn, false.positive = fp, false.negative = fn, 
				Accuracy = Accuracy, F1score = F1score, Precision = Precision, Recall = Recall, FPR = FPR)
	roc.matrix <- as.matrix(roc.list)
	rownames(roc.matrix) <- c("true positive", "true negative", "false positive", "false negative", 
				 "accuracy", "balanced F-score", "positive predictive value", "true positive rate", "false positive rate")
	return(roc.matrix)
}
# for comparing the result with different packages or approaches according to the true graph
compare = function (true.g, est.g, est.g2 = NULL, colnames = NULL) 
{
compare.true <- roc(true.g = true.g, est.g = true.g)
compare.g <- roc(true.g = true.g, est.g = est.g)
if (!is.null(est.g2)){
   compare.g2 <- roc(true.g = true.g, est.g = est.g2)
   compare.all <- cbind(compare.true, compare.g, compare.g2)
   if (is.null(colnames)){
     colnames <- c("True graph", "est.g", "est.g2")
     }
   colnames(compare.all) <- colnames
   return(compare.all)
   } else {
    compare.all <- cbind(compare.true, compare.g)
	if (is.null(colnames)){
     colnames <- c("True graph", "est.g")
     }
	colnames(compare.all) <- colnames
    return(compare.all)
   }
}
# Data generator according to the graph structure
bdgraph.sim = function(n = 1, p = 10, graph = "random", prob = NULL, v = NULL, u = NULL, A = NULL, 
                             K = NULL, sigma = NULL, vis = FALSE)
{
	if(graph == "random" & is.null(A)){
		if(is.null(prob)) prob <- 0.2
		A <- matrix(0, p, p)
		A[upper.tri(A)] <- rbinom(p * (p - 1) / 2, 1, prob)
		A <- A + t(A)
	}
	
	if (!is.null(A)) diag(A) <- 0
	
    if(is.null(sigma) & is.null(K)){
        if(is.null(u)) u <- 0.1
	    if(is.null(v)) v <- 0.3	  	
		K <- A * v
		# make K positive definite and standardized
		diag(K) <- abs(min(eigen(K) $ values)) + 0.1 + u
		sigma <- cov2cor(solve(K))
		K <- solve(sigma)
	} else {
	    if (is.null(sigma)) sigma <- solve(K)
	}
	
	# generate multivariate normal data
	d <- mvrnorm(n, rep(0, p), sigma)
	# graph and covariance visulization
	if (vis == TRUE){
	   g <- network(A, directed = F)
       plot.network(g, label = network.vertex.names(g), mode = "circle", main = "Graph structure")
	}
	return(list(data = d, sigma = sigma, K = K, A = Matrix(A, sparse = TRUE)))
}
# plot for class bdgraph
plot.bdgraph = function(x, g = 1, mode = "circle", edge.col = "black", label.col = "black", vertex.col = "red", ...)
{
  list.A <- x $ As
  lambda <- x $ lambda
  p <- x $ p
  prob.A <- lambda / sum(lambda)
  graphi <- list()
  gv <- c(rep(0, p * (p - 1) / 2))
  if (g == 2) op <- par(mfrow = c(1, 2), pty = "s")
  if (g > 2 & g < 7)  op <- par(mfrow = c(2, g %% 2 + trunc(g /2)), pty = "s")
  for (i in 1 : g){
    if (g > 6) dev.new()  
    gi <- list.A[[which(prob.A == sort(prob.A, decreasing = T)[i])]]
    gv <- 0 * gv
    gv[which(unlist(strsplit(as.character(gi), "")) == 1)] = 1
    graphi[[i]] <- matrix(0, p, p)
    graphi[[i]][upper.tri(graphi[[i]])] <- gv
    G <- network(graphi[[i]], directed = F)
    if (i == 1){
	   if (g == 1){
         main = paste(c("BEST GRAPH with size ", sum(graphi[[1]]), " : graph with highest probability"), collapse = "") 
        } else { 
		 main = paste(c("BEST GRAPH with size ", sum(graphi[[1]])), collapse = "")   
		   }
	} else {
	   main = paste(c(i, "th graph : size=", sum(graphi[[i]])), collapse = "")
	   }
    plot.network(G, label = network.vertex.names(G), mode = mode, main = main,
	 edge.col = edge.col, label.col = label.col, vertex.col = vertex.col,
     sub = paste(c("Posterior probability of graph = ", round(sort(prob.A, decreasing = TRUE)[i], 4)), collapse = ""))
  }
  if (g > 1 & g < 7)  par(op)
}
# summary of the result according to bdgraph
summary.bdgraph = function(object, plot = TRUE, mode = "circle", edge.col = "black", label.col = "black", vertex.col = "red", ...)
{
  As <- object $ As
  lambda <- object $ lambda
  p <- object $ p
  gv <- c(rep(0, p * (p - 1) / 2))
  graphi <- matrix(0, p, p)
  prob.A <- lambda / sum(lambda)
  gi <- As[[which(prob.A == max(prob.A))]]
  gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
  graphi[upper.tri(graphi)] <- gv
  if (plot == TRUE){
	# plot best graph
	G <- network(graphi, directed = F)
	if (!is.null(object $ allA)){
	   op <- par(mfrow = c(2, 2), pty = "s")
    } else {
        op <- par(mfrow = c(1, 2), pty = "s")
       }	
	plot.network(G, label = network.vertex.names(G), mode = mode, main = paste(c("BEST GRAPH with size ", sum(graphi)), collapse = ""),
	edge.col = edge.col, label.col = label.col, vertex.col = vertex.col,
	sub = paste(c("Posterior probability of graph = ", round(max(prob.A), 4)), collapse = ""))
	# plot posterior distribution of graph size
	suma <- sapply(As, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
	xx <- unique(suma)
	lambdag <- vector()
	for (i in 1:length(xx)){
	  lambdag[i] <- sum(lambda[which(suma == xx[i])])
	}
	plot(x = xx, y = lambdag, type = "h", main = "Posterior probability of graph size",
	ylab = "Pr(number of links in the graph|data)", xlab = "number of links in the graph")  
	if (!is.null(object $ allA)){
		# plot trace of graph size
		if (!is.null(object $ allA)){
		 allA <- object $ allA
		 y <- sapply(allA, function(x) length(which(unlist(strsplit(as.character(x), "")) == 1)))
		 plot(x = (object $ alla) * (1 : length(allA)), y, type = "l", main = "trace for size of graphs",
			  ylab = "sum of links in the graphs", xlab = "iterations")	
		}
		# plot ACF
		acf(y, main = "AIC for graph sizes") 
	}
	par(op)
  }
  # phat
  pvec <- 0 * gv
  for (i in 1:length(As)){
	inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
	pvec[inp] <- pvec[inp] + lambda[i]
	}
  phat <- 0 * graphi
  phat[upper.tri(phat)] <- pvec / sum(lambda)
  # estimation for precision matrix 
  Khat <- object $ Khat
  return.list <- list(best.graph = Matrix(graphi + t(graphi), sparse = TRUE), phat = Matrix(round(phat, 2), sparse = TRUE), 
  Khat = round(Khat, 3))
  return(return.list)
}
# print of the result according to bdgraph
print.bdgraph = function(x, round = 3, Khat = FALSE, phat = FALSE, ...)
{
  As <- x $ As
  lambda <- x $ lambda
  p <- x $ p
  # best graph
  prob.A <- lambda / sum(lambda)
  gv <- c(rep(0, p * (p - 1) / 2))
  gi <- As[[which(prob.A == max(prob.A))]]
  gv[which(unlist(strsplit(as.character(gi), "")) == 1)] <- 1
  graphi <- matrix(0, p, p)
  graphi[upper.tri(graphi)] <- gv
  
  cat(paste(""), fill = TRUE)
  cat(paste("Adjacency matrix of best graph"), fill = TRUE)
  cat(paste(""), fill = TRUE)
  printSpMatrix(Matrix(graphi + t(graphi), sparse = TRUE), col.names = TRUE, note.dropping.colnames = FALSE)
  
  cat(paste(""), fill = TRUE)
  cat(paste("Size of best graph =", sum(graphi)), fill = TRUE)
  cat(paste("Posterior probability of best graph = ", round(max(lambda) / sum(lambda), round)), fill = TRUE)  
  cat(paste(""), fill = TRUE)
  # print for precisiton matrix
  if (Khat == TRUE){
	  cat(paste(""), fill = TRUE)
	  cat(paste("Estimation of precision matrix"), fill = TRUE)
	  cat(paste(""), fill = TRUE)
	  print(round(x $ Khat, round))
	}
  # print for phat
  if (phat == TRUE){
	  pvec <- 0 * gv
	  for (i in 1:length(As)){
		 inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
		 pvec[inp] <- pvec[inp] + lambda[i]
		}
	  phat <- 0 * graphi
	  phat[upper.tri(phat)] <- pvec / sum(lambda)
	  cat(paste(""), fill = TRUE)
	  cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
	  cat(paste(""), fill = TRUE)
	  printSpMatrix(Matrix(round(phat, round), sparse = TRUE), col.names = TRUE, note.dropping.colnames = FALSE)  
    }
} 









