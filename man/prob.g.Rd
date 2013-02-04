\name{prob.g}
\alias{prob.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior probability for a specific graph
}
\description{
According to the output of the BDMCMC algorithm, this function gives us the 
posterior probability for a specific graph.
}
\usage{
prob.g(A, output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
upper triangular matrix in which \eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}.}
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
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
  truK # precision matrix of the true graph
  
  # generate the data (200 observations) from multivariate normal 
  # distribution with mean zero and percision matrix "truK"
  data <- mvrnorm(200, c(rep(0,p)), solve(truK))  
  output <- bdmcmc(data, meanzero = T, iter = 2000)
  
  A <- ceiling (truK)
  A[lower.tri[A, biag = T]] <- 0
  A # matrix A shows circle graph with 8 links and 8 nodes which is the true graphical model
  
  prob.g(A, output)
  }
}
