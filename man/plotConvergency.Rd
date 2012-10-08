\name{plotConvergency}
\alias{plotConvergency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative occupancy fractions for checking the convergency
}
\description{
Plot the cumulative occupancy fractions for all possible edges. 
It is a tool for checking the convergency of the BDMCMC algorithm.
}
\usage{
plotConvergency(output, skip = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{skip}{
it is a option for getting fast result for a cumulative plot according to part of the iterations.
}
}
\details{
Note that a spending time for this function depences on the graph. It should be slow for the 
high-dimensional graphs. To make it faster you can choose bigger value for 'skip', as you 
can see in following example.
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, arXiv:1210.5371. \url{http://arxiv.org/abs/1210.5371}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}

\seealso{
\code{\link{bdmcmc.low}}
}
\examples{
\dontrun{
  require(MASS)
  
  p <- 8 # number of nodes 
  # "truK" is the precision matrix of true graph
  truK <- diag(p)
  for (i in 1:(p-1)) truK[i,i+1] = truK[i+1,i] = 0.5
  truK[1,p] <- 0.4
  truK[p,1] <- 0.4
  
  # generate the data (200 observations) from multivariate normal 
  # distribution with mean zero and percision matrix "truK"
  data <- mvrnorm(200, c(rep(0,p)), solve(truK))  
  output <- bdmcmc(data, meanzero=T, iter=2000, print=T)
  
  # we run it for skip=5. For skip=1, it takes more time.
  plotConvergency(output, skip=5)
  }
}
