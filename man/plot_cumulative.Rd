\name{plot_cumulative}
\alias{plot_cumulative}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative occupancy fractions
}
\description{
Tools for checking the convergency of the BDMCMC algorithm
}
\usage{
plot_cumulative(output, distance = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
A list which is the result of BDMCMC algorithm 
}
  \item{distance}{
It is the option for fast result for a cumulative plot according to part of the iterations.
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
\code{\link{BDMCMC.low}}
}
\examples{
n=100; p=5
truK=diag(p)
for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
truK[1,p]=truK[p,1]=0.4
data=mvrnorm(n,c(rep(0,p)),solve(truK))
\dontshow{
  data(output.low)
  output <- output.low}
\dontrun{output = BDMCMC(data) }
# we run it for distance=20, you can run it for distance=1 and it`s less than 10 minites
plot_cumulative(output,distance=20)
}
