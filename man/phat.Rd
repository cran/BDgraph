\name{phat}
\alias{phat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior edge inclusion probabilities
}
\description{
According to the output of BDMCMC algorithm this function gives us the posterior edge inclusion 
probabioities for all possible edges.
}
\usage{
phat(output, round = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
A list which is the result of BDMCMC algorithm from the 'bdmcmc', 'bdmcmc.low', or 'bdmcmc.high' functions. 
}
  \item{round}{
A number for rounding all probabilities to the specified number of decimal places (default 3).
}
}
\value{
\item{phat}{Upper truculer matrix which shows the posterior inclusion probabilities for all possible links.}

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
  output = bdmcmc(data,meanzero=T,print=T)
  phat(output,2)
  }
}

