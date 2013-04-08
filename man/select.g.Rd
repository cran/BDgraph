\name{select.g}
\alias{select.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Selecting the best graphical models based on BDMCMC algorithm
}
\description{
According to output of the BDMCMC algorithm, this function gives us the
 graphical models with the highest posterior probabilities.  
}
\usage{
select.g(output, g = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of the BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{g}{
you can select the number of graphical models with highest probabilities (default is 1).
}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdmcmc}} and \code{\link{prob.allg}}
}
\examples{
\dontrun{
  p <- 8 # number of nodes 
  # "truK" is the precision matrix of true graph
  truK <- diag(p)
  for (i in 1:(p-1)) truK[i,i+1] <- truK[i+1,i] <- 0.5
  truK[1,p] <- truK[p,1] <- 0.4
  truK # precision matrix of the true graph
  
  # generate the data (200 observations) from multivariate normal 
  # distribution with mean zero and percision matrix "truK"
  data <- mvrnorm(200, c(rep(0,p)), solve(truK))
  
  output <- bdmcmc(data, meanzero = T, iter = 2000)
  select.g(output)
  }
}
