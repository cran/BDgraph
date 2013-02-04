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
prob.allg(output, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{verbose}{
logical: if verbose = TRUE, tracing information printing is abled. The default value is TRUE.
}
}
\value{
\item{list.A }{ a list which includes all the grpahs that the BDMCMC algorithm visits them.}
\item{prob.A }{ a vector which includes posterior probabilities of all graphs in 'list.A'.}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
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
  for (i in 1:(p-1)) truK[i,i+1] <- truK[i+1,i] <- 0.5
  truK[1,p] <- truK[p,1] <- 0.4
  truK
  
  # generate the data (200 observations) from multivariate normal 
  # distribution with mean zero and percision matrix "truK"
  data <- mvrnorm(200, c(rep(0,p)), solve(truK))  
  
  # First, we run the 'bdmcmc' function for small number of iterations
  output <- bdmcmc(data, iter = 40, burn = 30, meanzero = T)
  prob.allg(output)
  }
}
