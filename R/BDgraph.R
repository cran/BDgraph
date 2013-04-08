## All functions for the "BDgraph" package
###########################################################
## function for only obtaining elements of Psi matrix
Psi = function(A, b, H, p)
{
  nu <- apply(A, 1, sum)
  psi <- diag(sqrt(rchisq(p, b + nu)))
  psi[A == 1] <- rnorm(1)
  for (i in 1:(p-1)){
    for (j in (i+1):p){
        if (A[i,j] == 0){
            psi[i,j] <- - sum(psi[i, i:(j-1)] * H[i:(j-1), j])
            if (i > 1){
               for (r in 1:(i-1)){
                 psi[i,j] <- psi[i,j] - ((sum(psi[r, r:i] * H[r:i, i])) *
                             (sum(psi[r, r:j] * H[r:j, j]))) / (psi[i,i])
               }
            }
        }
    }
  }
  return(psi)
}
# Algorithm 3.1: function for Monte Carlo approxiation for expection in normalizing constant
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
# function for computing Normalizing constans of G-Wishart distribution according to ATAY-KAYIS AND mASSAM (2005)
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
   for (i in 1:p){
      c_dT <- c_dT + ((sum(A[i,]) / 2) * log(pi) + ((b + 2 * sum(A[i,])) / 2) * log(2) +
      lgamma((b + sum(A[i,])) / 2) + (b + sum(A[i,]) + sum(A[,i])) * log(Ti[i,i]))
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
   cont <- 1
   for (i in 1:(iter * 1000)){  
      psi.new <- Psi(A, b, H, p)
	  alpha <- exp((sum((1 - (diag(p) + A)) * psi * psi) - 
	  sum((1 - (diag(p) + A)) * psi.new * psi.new)) / 2)
	  if (runif(1) < alpha){ 
	    cont <- cont + 1
		psi.new <- psi
	    }  
      if (cont == iter) break
	}
    return(t(psi %*% Ts) %*% (psi %*% Ts))
}
# Sampling for G-Wishart distribution according to Blocked Gibbs sampling, Wang (2012)
block.gibbs = function (K, A, bstar, Ds, gibbs.iter = 1)
{
p <- ncol(A)
for (k in 1 : gibbs.iter){
  for (i in 1 : (p - 1)){
	if (sum(A[i,]) != 0){
      for (j in (i + 1) : p){
	    if (A[i,j] == 1){
		  B <- Ds[c(i,j),c(i,j)]
		  a <- rWishart(1, df = bstar + (p - 2), Sigma = solve(B))
		  a <- matrix(a, 2, 2)
		  Kvv <- K[c((1 : p)[- c(i,j)]), c((1 : p)[- c(i,j)])]
          kjv <- K[c(i,j),c((1:p)[- c(i,j)])]
          Kc <- a + (kjv) %*% (solve(Kvv)) %*% t(kjv)
          K[c(i,j),c(i,j)] <- (Kc + t(Kc)) / 2 # Numerical stable
		}
      }
	}
    if (sum(A[i,]) + sum(A[,i]) == 0){
	  a <- rgamma(1, (bstar + (p - 1)) / 2, Ds[i,i] / 2)
	  Kvv <- K[c((1 : p)[- i]), c((1 : p)[- i])]
      kjv <- matrix(K[i, c((1 : p)[- i])], 1, (p - 1))
      K[i,i] <- a + (kjv) %*% (solve(Kvv)) %*% t(kjv)
	}
  }
  if (sum(A[,p]) == 0){
	a <- rgamma(1, (bstar + (p - 1)) / 2, Ds[p,p] / 2)
	Kvv <- K[c((1 : p)[- p]), c((1 : p)[- p])]
    kjv <- matrix(K[p,c((1:p)[-p])], 1, (p - 1))
    K[p,p] <- a + (kjv) %*% (solve(Kvv)) %*% t(kjv)
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
bdmcmc.low = function(data, n = NULL, meanzero = FALSE, iter = 5000, burn = floor(iter / 2), skip = 1, 
gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, A = "full", MCiter = 10, summary = FALSE, verbose = TRUE, all.A = FALSE)
{
  if (is.matrix(data) == F & is.data.frame(data) == F){
     stop("Data should be a matrix or dataframe")
  }
  if (is.data.frame(data) == T) data <- data.matrix(data)
  if (any(is.na(data))) stop("Data should contain no missing data") 
  if (iter <= burn){
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
  start.A <- pmatch(A, c("full", "empty", "glasso"))[1]	
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
  if (is.null(A)) {
     A <- 0 * S
     A[upper.tri(A)] <- 1
     } 
  if (is.null(D)) D <- diag(p)
  bstar <- b + n
  Ds <- D + S
  Dsinv <- solve(Ds)
  Ts <- chol(Dsinv)
  Ti <- chol(D)
  H <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
  Hs <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))
  K <- sampleK(A, bstar, Hs, Ts, p)
  As <- allA  <- vector() # vector of numbers like "10100"
  lambda <- all.lambda <- vector() # waiting time for every state
  alla <- ceiling(iter / 2000)# for saving all.A which we need it for plotConvergency function
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
	if ( all.A == TRUE & g %% alla == 0){
        indA <- paste(A[upper.tri(A)], collapse = '')
		allA <- c(allA, indA)
		all.lambda <- c(all.lambda, 1 / sum(rates))
    }
    if (g > burn && g %% skip == 0){
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
    K <- block.gibbs(K = K, A = A, bstar = bstar, Ds = Ds, gibbs.iter = 1)
  }
  if(verbose == TRUE){
	mes = paste(c("    ", iter," iterations done.                              "), collapse = "")
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }  
  if (summary == FALSE){
     if (all.A == TRUE){
        return(list(As = As, lambda = lambda, p = p, allA = allA, all.lambda = all.lambda, alla = alla))
	 } else {
	    return(list(As = As, lambda = lambda, p = p))
	 }
  } else {
        output <- list(As = As, lambda = lambda, p = p)
		# best graph and estimaition of its parameters
		bestg <- select.g (output, g = 1)
		print(round(bestg, 2))
        # for phat
        pvec <- c(rep(0, p * (p - 1) / 2))
        for (i in 1:length(As)){
           inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
	       pvec[inp] <- pvec[inp] + lambda[i]
        }
        phat <- matrix(0, p, p)
        phat[upper.tri(phat)] <- pvec / sum(lambda)
		cat(paste(""), fill = TRUE)
		cat(paste(""), fill = TRUE)
		cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
		cat(paste(""), fill = TRUE)
        return(round(phat, 2))	
    }
}
# function for obtaining "c" (Wang2012page182) from matrix K
get.c = function(K, p, j){
   kjj <- K[c((1:p)[- j]), c((1:p)[- j])]
   kjv <- matrix(K[j,][- j], 1, (p - 1))
   return((kjv) %*% (solve(kjj)) %*% (t(kjv)))  
}
# function for computing "c.star - c" from the matrix K for k_jj
get.cs_c = function(K, p, i, j){
   kjj <- K[c((1:p)[- j]), c((1:p)[- j])]
   kjv <- K[j,][- j]
   B <- solve(kjj)
   return(K[i,j] * (K[i,j] * B[i,i] - 2 * (kjv %*% B[,i])))  
}
## Algorithm 2.1: BD-MCMC algorithm for high-dimentional problem (roughly graphs with more than 8 nodes)
bdmcmc.high = function(data, n = NULL, meanzero = FALSE, iter = 5000, burn = floor(iter / 2),
skip = 1, gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, A = "full", summary = FALSE, verbose = TRUE, all.A = FALSE)
{
  if (is.matrix(data) == F & is.data.frame(data) == F){
     stop("Data should be a matrix or dataframe")
  }
  if (is.data.frame(data) == T) data <- data.matrix(data)
  if (any(is.na(data))) stop("Data should contain no missing data") 
  if (iter <= burn){
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
  start.A <- pmatch(A, c("full", "empty", "glasso"))[1]	
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
  if (is.null(D)) D <- diag(p)
  Ti <- chol(solve(D))
  bstar <- b + n
  Ds <- D + S
  Ts <- chol(solve(Ds))
  H <- Ti / t(matrix(rep(diag(Ti) ,p), p, p))
  Hs <- Ts / t(matrix(rep(diag(Ts) ,p), p, p))
  K <- sampleK(A, bstar, Hs, Ts, p, iter = 10)
  As <- allA  <- vector() # vector of numbers like "10100"
  lambda <- all.lambda <- vector() # waiting time for every state
  alla <- ceiling(iter / 2000) # for saving allA which we need it for plotConvergency function
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
	if ( all.A == TRUE & g %% alla == 0){
        indA <- paste(A[upper.tri(A)], collapse = '')
		allA <- c(allA, indA)
		all.lambda <- c(all.lambda, 1 / sum(rates))
    }
    if (g > burn && g %% skip == 0){
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
	K <- block.gibbs(K, A, bstar, Ds, gibbs.iter = 1)
  }
  if(verbose == TRUE){
	mes = paste(c("    ", iter," iterations done.                              "), collapse = "")
    cat(mes, "\r")
    cat("\n")
    flush.console()
  }  
  if (summary == FALSE){
     if (all.A == TRUE){
        return(list(As = As, lambda = lambda, p = p, allA = allA, all.lambda = all.lambda, alla = alla))
	 } else {
	    return(list(As = As, lambda = lambda, p = p))
	 }
  } else {
        output <- list(As = As, lambda = lambda, p = p)
		# best graph and estimaition of its parameters
		bestg <- select.g (output, g = 1)
		print(round(bestg, 2))
        # for phat
        pvec <- c(rep(0, p * (p - 1) / 2))
        for (i in 1:length(As)){
           inp <- which(unlist(strsplit(as.character(As[i]), "")) == 1)
	       pvec[inp] <- pvec[inp] + lambda[i]
        }
        phat <- matrix(0, p, p)
        phat[upper.tri(phat)] <- pvec / sum(lambda)
		cat(paste(""), fill = TRUE)
		cat(paste(""), fill = TRUE)
		cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
		cat(paste(""), fill = TRUE)
        return(round(phat, 2))	
    }
}
## Main function: BDMCMC algorithm for selecting the best graphical model
bdmcmc = function(data, n = NULL, meanzero = FALSE, iter = 5000, burn = floor(iter / 2), skip = 1,
gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, A = "full", MCiter = 10, summary = FALSE, verbose = TRUE, all.A = FALSE)
{
  p <- ncol(data)
  if (p < 8){
    return(bdmcmc.low(data, n, meanzero, iter, burn, skip, gamma.b, prior.g, b, D, A, MCiter, summary, verbose, all.A))
  } else {
    return(bdmcmc.high(data, n, meanzero, iter, burn, skip, gamma.b, prior.g, b, D, A, summary, verbose, all.A))
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
# function for computing the probability of one especial graph
prob.g = function(A, output)
{
   As <- output $ As
   lambda <- output $ lambda
   indA <- paste(A[upper.tri(A)], collapse = '')
   lambda.g <- 0
   for (i in 1:length(As)){
       if (identical(indA,As[i]) == TRUE) lambda.g <- lambda[i]
	   }
   cat(paste(""), fill = TRUE)
   mes = paste(c("  Posterior probability of the graph = ", round(lambda.g / sum(lambda), 4)), collapse = "")
   cat(mes, "\r")
   cat("\n")
}
# plot for probability of graphs according to number of their links
plotLinks = function(output, xlim = c(0, output $ p * (output $ p - 1) / 2))
{
   p <- output $ p
   As <- output $ As
   lambda <- output $ lambda
   nominator <- c(rep(0, xlim[2] - xlim[1] + 1))
   for (i in 1:length(As)){
      for (j in xlim[1]:xlim[2]){
	      lAs <- length(which(unlist(strsplit(as.character(As[i]), "")) == 1))
          if (lAs == j) {nominator[j - xlim[1] + 1] <- nominator[j - xlim[1] + 1] + lambda[i]}
      }
   }
   plot(x = xlim[1]:xlim[2], y = nominator / sum(lambda), type = "h", main = "",
   ylab = "Pr(number of links in the graph|data)", xlab = "number of links in the graph")
}
# function for checking the convergency of the BDMCMC algorithm
plotConvergency = function(output, skip = ceiling(length(output $ allA) / 2000), verbose = TRUE)
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
  ffv <- 0 * ff[,1]
  ffv[inp] <- all.lambda.new[1]
  for (g in 2:length.allA.new){
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
  abline(v = skip * (output $ alla) * length.allA.new / 2, col = "blue")
}
# plot sum of the links in the graphs for checking the convergency of BDMCMC algorithm
plotSum = function(output, xlim = c(0, length(output $ allA)), main = "")
{
    allA <- output $ allA
    iter <- length(allA)
    y <- 0 * (xlim[1]:xlim[2])
    for (i in (xlim[1] + 1):xlim[2]){
         y[i - xlim[1] + 1] <- length(which(unlist(strsplit(as.character(allA[i]), "")) == 1))
    }
    plot(x = (output $ alla) * (xlim[1]:xlim[2]), y, type = "l",
    ylab = "sum of links in the graphs", xlab = "iterations")
}
# Program for computing the probability of all possible graphical models by using the result of BDMCMC algorithm
prob.allg = function(output)
{
  list.A <- output $ As
  lambda <- output $ lambda
  return(list(list.A = list.A, prob.A = lambda / sum(lambda)))
}                                              
# function summary of the result
select.g = function (output, g = 1)
{
  list.A <- output $ As
  lambda <- output $ lambda
  p <- output $ p
  prob.A <- lambda / sum(lambda)
  best.graph <- list.A[[which(prob.A == max(prob.A))]]
  # transferring to value likes "10100" to matrix A
  gv <- c(rep(0, p * (p - 1) / 2))
  gv[which(unlist(strsplit(as.character(best.graph), "")) == 1)] = 1
  best.graph <- matrix(0, p, p)
  best.graph[upper.tri(best.graph)] <- gv
  for (i in 1:g){
    dev.new()
    graphi <- list.A[[which(prob.A == sort(prob.A, decreasing = T)[i])]]
    gv <- c(rep(0, p * (p - 1) / 2))
    gv[which(unlist(strsplit(as.character(graphi), "")) == 1)] = 1
    graphi <- matrix(0, p, p)
    graphi[upper.tri(graphi)] <- gv
    G <- network(graphi, directed = F)
    if (i == 1){
         main = "BEST GRAPH: graph with highest probability"
    } else {
      if (i == 2){
        main = "Graph with 2nd highest probability"
      } else {
        if (i == 3){
          main = "Graph with 3rd highest probability"
        } else {
          main = paste(c("Graph with ", i, "th highest probability"), collapse = "")
        }
      }
    }
    plot.network(G, label = network.vertex.names(G), main = main,
    sub = paste(c("Posterior probability of graph = ", round(sort(prob.A, decreasing = TRUE)[i], 3)), collapse = ""))
  }
  cat(paste(""), fill = TRUE)
  cat(paste("Adjacency matrix of best graph"), fill = TRUE)
  cat(paste(""), fill = TRUE)
  return(best.graph)
}
# for comparing the result with different packages or approaches
compare = function (est.g, true.g) 
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
compare.list <- list(true.positive = tp, true.negative = tn, false.positive = fp, false.negative = fn, 
            Accuracy = Accuracy, F1score = F1score, Precision = Precision, Recall = Recall, FPR = FPR)
compare.matrix <- as.matrix(compare.list)
rownames(compare.matrix) <- c("true positive", "true negative", "false positive", "false negative", 
             "accuracy", "balanced F-score", "positive predictive value", "true positive rate", "false positive rate")
return(compare.matrix)
}







