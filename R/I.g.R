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
