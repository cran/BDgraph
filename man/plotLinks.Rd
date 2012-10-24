\name{plotLinks}
\alias{plotLinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of posterior distribution for graphs according to number of their links
}
\description{
This function is for ploting the posterior distribution of all possible graphs according to number
of their links.
}
\usage{
plotLinks(output, xlim = c(0, (nrow(output$As[[1]])) * (nrow(output$As[[1]]) - 1)/2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{xlim}{
the x limits (min, max) of the plot.
}
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, arXiv:1210.5371. \url{http://arxiv.org/abs/1210.5371}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}

\seealso{
\code{\link{bdmcmc}}
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
  
  output <- bdmcmc(data, meanzero=T, iter=500, print=T)
  plotLinks(output)
  }
}
