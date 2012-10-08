\name{prob.g}
\alias{prob.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior probability for one special graph
}
\description{
According to the output of BDMCMC algorithm, this function gives us the 
posterior probability for one special graphical model.
}
\usage{
prob.g(A, output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Upper truculer matrix in which a_ij=1 if there is a link between notes i and j, otherwise a_ij=0.
} 
  \item{output}{
A list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
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
\code{\link{bdmcmc}}
}
\examples{
\dontrun{
  require(MASS)
  n=200; p=8
  truK=diag(p)
  for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
  truK[1,p]=truK[p,1]=0.4
  data=mvrnorm(n,c(rep(0,p)),solve(truK))
  output = bdmcmc(data,meanzero=T) 
  A=0*truK
  for (i in 1:(p-1)){
    for (j in (i+1):p){
       if (truK[i,j]!=0) {A[i,j]=1}
    }
  }
  A[1,p]=1
  prob.g(A, output)
  }
}
