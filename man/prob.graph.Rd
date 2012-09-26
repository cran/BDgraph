\name{prob.graph}
\alias{prob.graph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior probability of the graph
}
\description{
According to the output of BDMCMC algorithm, this function gives us the 
posterior probability of the graphical models.
}
\usage{
prob.graph(A, output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Upper truculer matrix which shows the starting graphs in which a_ij=1 if there 
is a link between notes i and j, otherwise a_ij=0  
}
  \item{output}{
  A list which is the result of BDMCMC algorithm
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
\code{\link{BDMCMC}}
}
\examples{
n=200; p=8
truK=diag(p)
for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
truK[1,p]=truK[p,1]=0.4
data=mvrnorm(n,c(rep(0,p)),solve(truK))
\dontshow{
  data(output.high)
  output <- output.high}
\dontrun{output = BDMCMC(data) }
A=0*truK
for (i in 1:(p-1)){
    for (j in (i+1):p){
       if (truK[i,j]!=0) {A[i,j]=1}
       }
    }
A[1,p]=1
prob.graph(A, output)
}
