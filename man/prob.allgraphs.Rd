\name{prob.allgraphs}
\alias{prob.allgraphs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior probability of all possible grpahs
}
\description{
According to the output of BDMCMC algorithm this function gives us the 
posterior probability of all possible graphical models.
}
\usage{
prob.allgraphs(output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
  A list which is the result of BDMCMC algorithm
}
}
\value{
\item{list.A }{ A list which include all the grpahs that BDMCMC algorithm visit them}
\item{prob.A }{ A vector which include posterior probabilities of graphs in list.A}
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
\code{\link{BDMCMC}}
}
\examples{
n=200; p=8
truK=diag(p)
for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
truK[1,p]=truK[p,1]=0.4
data=mvrnorm(n,c(rep(0,p)),solve(truK))
# we run our BDMCMC for small number of iterations, only for checking the result
output=BDMCMC(data,iter=40,burn=10)
prob.allgraphs(output)
}
