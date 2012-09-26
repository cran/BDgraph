\name{Phat}
\alias{Phat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior edges inclusion probabilities
}
\description{
According to the output of BDMCMC algorithm this function gives us the posterior 
probabioity of all possible links in the graph.
}
\usage{
Phat(output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
A list which is the result of BDMCMC algorithm 
}
}
\value{
\item{Phat}{Upper truculer matrix which shows the posterior inclusion probabilities for all
the possible links in the graph}

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
round(Phat(output),2)
}

