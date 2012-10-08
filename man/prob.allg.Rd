\name{prob.allg}
\alias{prob.allg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior probability of all possible graphs
}
\description{
According to the output of the BDMCMC algorithm, this function gives us the 
posterior probability of all possible graphical models. Aslo, it give us
all graphs that the BDMCMC algorithm visits them.
}
\usage{
prob.allg(output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
}
\value{
\item{list.A }{ a list which includes all the grpahs that the BDMCMC algorithm visits them.}
\item{prob.A }{ a vector which includes posterior probabilities of all graphs in 'list.A'.}
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
  
  # First, we run the 'bdmcmc' function for small number of iterations
  output <- bdmcmc(data, iter=40, burn=30, meanzero=T, print=T)
  prob.allg(output)
  }
}
