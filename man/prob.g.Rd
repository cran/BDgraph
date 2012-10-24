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
  output <- bdmcmc(data, meanzero=T, iter=1000, print=T)
  
  A <- 0*truK
  for (i in 1:(p-1)){
    for (j in (i+1):p){
       if (truK[i,j]!=0) A[i,j] <- 1
    }
  }
  A[1,p] <- 1
  # matrix A shows circle graph with 8 links and 8 nodes which is the true graphical model
  prob.g(A, output)
  }
}
