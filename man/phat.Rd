\name{phat}
\alias{phat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior edge inclusion probabilities
}
\description{
According to the output of BDMCMC algorithm, this function gives us the posterior edge inclusion 
probabilities for all possible edges.
}
\usage{
phat(output, round = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{round}{
a number for rounding all probabilities to the specified number of decimal places (default 3).
}
}
\value{
\item{phat}{upper triangular matrix which shows the posterior inclusion probabilities for all possible edges.}

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
  phat(output, round=2)
  }
}

