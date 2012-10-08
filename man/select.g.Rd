\name{select.g}
\alias{select.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Selecting the best graphical models based on BDMCMC algorithm
}
\description{
According to output of BDMCMC algorithm, this function gives us the
 graphical models with the highest posterior probabioities.  
}
\usage{
select.g(output, g = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
A list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{g}{
You can select the number of graphical models with highest probabilities (default 2).
}
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, Journal of the Royal Statistical Society: Series B,
submitted.
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdmcmc}} and \code{\link{prob.allg}}
}
\examples{
\dontrun{
  require(MASS)
  n=200; p=8
  truK=diag(p)
  for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
  truK[1,p]=truK[p,1]=0.4
  data=mvrnorm(n,c(rep(0,p)),solve(truK))
  output = bdmcmc(data,meanzero=T,print=T)
  select.g(output)
  }
}
