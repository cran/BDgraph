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
	   dd[lower.tri(dd == 0, diag = T)] = 1
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
          kjv <- K[c(i,j),c((1:p)[-c(i,j)])]
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
gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, A = "full", MCiter = 10, summary = FALSE, verbose = TRUE)
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
  if (is.null(Sn$n) & is.null(n)){
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
		A <- as.matrix(Matrix(A$refit, sparse = F))
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
  Ks <- As <- rates_sample <- allA <- list()
  lambda <- vector() ## waiting time in every state
  cont <- allAcont <- 0
  alla <- ceiling(iter / 2000)# for saving allA which we need it for plotConvergency function
  for (g in 1:iter){
    if(verbose == T){
	   mes <- paste(c("    MCMC iterations : ", g, " from ", iter), collapse="")
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
        (sum(diag(Ds %*% (Kmi - K)))) / 2)
        if (is.finite(rates[i,j]) == FALSE) rates[i,j] <- gamma(170)
        }
     }
    }
    if (g %% alla == 0){
        allAcont <- allAcont + 1
        allA[[allAcont]] <- A
    }
    if (g > burn && g %% skip == 0){
        cont <- cont + 1
        Ks[[cont]] <- K
        As[[cont]] <- A
        lambda[cont] <- sum(rates)
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
     return(list(Ks = Ks, As = As, lambda = lambda, allA = allA, alla = alla))
  } else {
        output <- list(Ks = Ks, As = As, lambda = lambda, allA = allA, alla = alla)
		# best graph and estimaition of its parameters
		bestg <- select.g (output, K = TRUE, g = 1)
		print(round(bestg, 2))
        # for phat
		phat <- 0 * As[[1]]
        for (i in 1:(p-1)){
            for (j in (i+1):p){
                for (k in 1:length(As)){
                    phat[i,j] <- phat[i,j] + As[[k]][i,j] / lambda[k]
                }
                phat[i,j] <- phat[i,j] / (sum(1 / lambda))
            }
        }
		cat(paste(""), fill = TRUE)
		cat(paste(""), fill = TRUE)
		cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
		cat(paste(""), fill = TRUE)
        return(round(phat, 2))	
    }
}
# function for obtaining "c" (Wang2012page182) from matrix K
get.c = function(K, p, j){
   kjj <- K[c((1:p)[-j]),c((1:p)[-j])]
   kjv <- matrix(K[j,][- j], 1, (p - 1))
   return((kjv) %*% (solve(kjj)) %*% (t(kjv)))  
}
# function for computing "c.star - c" from the matrix K for k_jj
get.cs_c = function(K, p, i, j){
   kjj <- K[c((1:p)[-j]),c((1:p)[-j])]
   kjv <- K[j,][- j]
   B <- solve(kjj)
   return(K[i,j] * (K[i,j] * B[i,i] - 2 * (kjv %*% B[,i])))  
}
## Algorithm 2.1: BD-MCMC algorithm for high-dimentional problem (roughly graphs with more than 8 nodes)
bdmcmc.high = function(data, n = NULL, meanzero = FALSE, iter = 5000, burn = floor(iter / 2),
skip = 1, gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, A = "full", summary = FALSE, verbose = TRUE)
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
		A <- as.matrix(Matrix(A$refit, sparse = F))
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
  Ks <- As <- allA <- list()
  lambda <- vector() # waiting time for every state
  cont <- allAcont <- 0
  alla <- ceiling(iter / 2000)# for saving allA which we need it for plotConvergency function
  for (g in 1:iter){
    if(verbose == T){
	   mes <- paste(c("    MCMC iterations : ", g, " from ", iter), collapse="")
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
                     (sum(diag(Ds %*% (Kmi - K)))) / 2 )
          if (is.finite(rates[i,j]) == FALSE) rates[i,j] <- gamma(170)
        }
      }
    }
	if (g %% alla == 0){
        allAcont <- allAcont + 1
        allA[[allAcont]] = A
    }
    if (g > burn && g %% skip == 0){
        cont <- cont + 1
        Ks[[cont]] <- K
        As[[cont]] <- A
        lambda[cont] <- sum(rates)
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
     return(list(Ks = Ks, As = As, lambda = lambda, allA = allA, alla = alla))
  } else {
        output <- list(Ks = Ks, As = As, lambda = lambda, allA = allA, alla = alla)
		# best graph and estimaition of its parameters
		bestg <- select.g (output, K = TRUE, g = 1)
		print(round(bestg, 2))
        # for phat
		phat <- 0 * As[[1]]
        for (i in 1:(p-1)){
            for (j in (i+1):p){
                for (k in 1:length(As)){
                    phat[i,j] <- phat[i,j] + As[[k]][i,j] / lambda[k]
                }
                phat[i,j] <- phat[i,j] / (sum(1 / lambda))
            }
        }
		cat(paste(""), fill = TRUE)
		cat(paste(""), fill = TRUE)
		cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
		cat(paste(""), fill = TRUE)
        return(round(phat, 2))	
    }
}
## Main function: BDMCMC algorithm for selecting the best graphical model
bdmcmc = function(data, n = NULL, meanzero = FALSE, iter = 5000, burn = floor(iter / 2), skip = 1,
gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, A = "full", MCiter = 10, summary = FALSE, verbose = TRUE)
{
  p <- ncol(data)
  if (p < 8){
    return(bdmcmc.low(data, n, meanzero, iter, burn, skip, gamma.b, prior.g, b, D, A, MCiter, summary, verbose))
  } else {
    return(bdmcmc.high(data, n, meanzero, iter, burn, skip, gamma.b, prior.g, b, D, A, summary, verbose))
  }
}
# function for comuting probability of all links in graph
phat = function(output, round = 3)
{
   As <- output$As
   lambda <- output$lambda
   p <- nrow(As[[1]])
   phat <- 0 * As[[1]]
   for (i in 1:(p-1)){
      for (j in (i+1):p){
         for (k in 1:length(As)){
             phat[i,j] <- phat[i,j] + As[[k]][i,j] / lambda[k]
         }
         phat[i,j] <- phat[i,j] / (sum(1 / lambda))
      }
   }
   cat(paste(""), fill = TRUE)
   cat(paste("Posterior edge inclusion probabilities for all possible edges equal with"), fill = TRUE)
   cat(paste(""), fill = TRUE)
   return(round(phat, round))
}
# function for computing the probability of one especial graph
prob.g = function(A, output)
{
   As <- output$As
   lambda <- output$lambda
   lambda.g <- vector()
   for (i in 1:length(As)){
       if (identical(A,As[[i]]) == TRUE) lambda.g <- c(lambda.g, lambda[i])
   }
   round(sum(1 / lambda.g) / sum(1 / lambda), 4)
}
# plot for probability of graphs according to number of their links
plotLinks = function(output, xlim = c(0, (nrow(output$As[[1]])) * (nrow(output$As[[1]]) - 1) / 2))
{
   p <- nrow(output$As[[1]])
   As <- output$As
   lambda <- output$lambda
   nominator <- c(rep(0, xlim[2] - xlim[1] + 1))
   for (i in 1:length(As)){
      for (j in xlim[1]:xlim[2]){
          if (sum(As[[i]]) == j) {nominator[j - xlim[1] + 1] <- nominator[j - xlim[1] + 1] + 1 / lambda[i]}
      }
   }
   plot(x = xlim[1]:xlim[2], y = nominator / sum(1 / lambda), type = "h", main = "",
   ylab = "Pr(number of links in the graph|data)", xlab = "number of links in the graph")
}
# function for checking the convergency of the BDMCMC algorithm
plotConvergency = function(output, skip = 1, verbose = TRUE)
{
  allA <- output$allA
  p <- nrow(allA[[1]])
  allA.new <- list()
  for (i in 1:length(allA)){
     g <- i * skip
     if (g > length(allA)) break
     allA.new[[i]] <- allA[[g]]
  }   
  length.allA <- length(allA.new)
  ff <- matrix(0, p * (p - 1) / 2, length.allA)
  for (g in 1:length.allA){
     if (verbose == TRUE){
	    mes <- paste(c("Calculating cumulative occupancy fractions....in progress : ", floor(100 * g / length.allA), "%"), collapse = "")
	    cat(mes, "\r")
	    flush.console()	
       }  
     for (k in 1:g){
        con <- 0
        for (i in 1:(p-1)){
           for (j in (i+1):p){
              con <- con + 1
              ff[con,g] <- ff[con,g] + allA.new[[k]][i,j] / g
           }
        }
     }
  }
  if(verbose == TRUE){
	mes = paste(c("Calculating cumulative occupancy fractions....done.                   "), collapse = "")
    cat(mes, "\r")
    cat("\n")
    flush.console()
  } 
  matplot(x = skip * ((output$alla) * (1:length.allA)), y = t(ff), type = "l", lty = 1, col = 1,
  xlab = "number of iterations", ylab = "cumulative occupancy fractions for each links")
}
# plot sum of the links in the graphs for checking the convergency of BDMCMC algorithm
plotSum = function(output, xlim = c(0, length(output$allA)))
{
    allA <- output$allA
    iter <- length(allA)
    y <- 0 * (xlim[1]:xlim[2])
    for (i in (xlim[1] + 1):xlim[2]){
         y[i - xlim[1] + 1] <- sum(allA[[i]])
    }
    plot(x = (output$alla) * (xlim[1]:xlim[2]), y, type = "l", main = "",
    ylab = "sum of links in the graphs", xlab = "iterations")
}
# Program for computing the probability of all possible graphical models by using the result of BDMCMC algorithm
prob.allg = function(output, verbose = TRUE)
{
   As.lim <- output$As
   lambda <- output$lambda
   lambda.lim <- lambda
   list.A <- list.lambda <- list()
   fre <- vector()
   i <- 0
   while (length(As.lim) > 1){
      i <- i + 1
      if (verbose == TRUE){
	     mes <- paste(c("    Model selection .... in progress : ", floor(100 * (1 - length(As.lim) / length(output$As))), "%"), collapse = "")
	     cat(mes, "\r")
	     flush.console()	
        }  
      list.A[[i]] <- As.lim[[1]]
      list.lambda[[i]] <- vector()
      free <- 0
      As.lim2 <- As.lim
      lambda.lim2 <- lambda.lim
      for (g in 1:length(As.lim)){
         if (identical(list.A[[i]], As.lim[[g]]) == TRUE){
            list.lambda[[i]] = c(list.lambda[[i]], lambda.lim[g])
            As.lim2 <- As.lim2[- (g - free)]
		        lambda.lim2 <- lambda.lim2[- (g - free)]
		        free <- free + 1
         }
      }
      As.lim <- As.lim2
      lambda.lim <- lambda.lim2
      fre[i] <- free
   }
   prob.A <- vector()
   for (i in 1:length(list.A)){
       prob.A[i] <- sum(1 / list.lambda[[i]]) / sum(1 / lambda)
   }
   if(verbose == TRUE){
	 mes = paste(c("    Model selection .... done.                                  "), collapse = "")
     cat(mes, "\r")
     cat("\n")
     flush.console()
   }  
   return(list(list.A = list.A, prob.A = prob.A))
}                                              
# function summary of the result
select.g = function (output, g = 1, K = FALSE, verbose = TRUE)
{
  output.allg <- prob.allg(output, verbose = verbose)
  list.A <- output.allg$list.A
  prob.A <- output.allg$prob.A
  best.graph <- list.A[[which(prob.A == max(prob.A))]]
  for (i in 1:g){
    dev.new()
    G <- network(list.A[[which(prob.A == sort(prob.A, decreasing = T)[i])]], directed = F)
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
    sub = paste(c("Posterior probability of graph=", round(sort(prob.A, decreasing = TRUE)[i], 3)), collapse = ""))
  }
  cat(paste(""), fill = TRUE)
  cat(paste("Adjacency matrix of best graph"), fill = TRUE)
  cat(paste(""), fill = TRUE)
  return(best.graph)
  if (K == TRUE){
	 Ks <- output$Ks
	 As <- output$As
	 bestG <- list.A[[which(prob.A == max(prob.A))]]
	 bestK <- 0 * Ks[[1]]
	 for (g in 1:length(Ks)){
	    if (identical(bestG, As[[g]]) == TRUE){
		     bestK <- bestK + Ks[[g]]
	        }
	    }
	 cat(paste(""), fill = TRUE)
	 cat(paste(c("The best graph has ", sum(bestG), " edges."),collapse = ""), fill = TRUE)	
	 cat(paste(""), fill = TRUE)
	 cat(paste("estimation of precision matrix for the best graph:"), fill = TRUE)
	 cat(paste(""), fill = TRUE)
	 return(round(bestK / length(bestK), 3))
	}
}






