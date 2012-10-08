\name{plotConvergency}
\alias{plotConvergency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cumulative occupancy fractions for checking the convergency
}
\description{
Plot the cumulative occupancy fractions for all possible edges. 
It is a tool for checking the convergency of the BDMCMC algorithm.
}
\usage{
plotConvergency(output, distance = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
A list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{distance}{
It is the option for fast result for a cumulative plot according to part of the iterations.
}
}
\details{
Not that the spending time for this function depences on the number of nodes in the graph. 
As a result, it should be slow for the high-dimensional graphs. To make it faster you can choose bigger 
value for dentance, as you can see in following example.
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
\code{\link{bdmcmc.low}}
}
\examples{
\dontrun{
  require(MASS)
  n=100; p=5
  truK=diag(p)
  for (i in 1:(p-1)) truK[i,i+1]=truK[i+1,i]=0.5
  truK[1,p]=truK[p,1]=0.4
  data=mvrnorm(n,c(rep(0,p)),solve(truK))
  output = bdmcmc(data,meanzero=T)
  # we run it for distance=20. For distance=1, it takes around 10 minites
  plotConvergency(output,distance=20)
  }
}
