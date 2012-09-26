\name{BDgraph-package}
\alias{BDgraph-package}
\alias{BDgraph}
\docType{package}
\title{
Bayesian model selection in Gaussian graphical models based on BDMCMC algorithm
}
\description{
BDgraph package is a statistical tool for Bayesian model determination in 
Gaussian graphical models based on birth-death MCMC methodology.      
}
\details{
\tabular{ll}{
Package: \tab BDgraph\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-09-23\cr
License: \tab GPL(>=3)\cr
}
This package includes several functions. The main function is BDMCMC which is 
brith-death MCMC algorithm for Bayesian model selection in Gaussian graphical models.
}
\author{
Abdolreza Mohammadi <a.mohammadi@rug.nl> and Ernst Wit
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, Journal of the Royal Statistical Society: Series B,
submitted.

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable gaussian graphical models. Biometrika 
92(2), 317-335.
}
\keyword{ package }
\seealso{
\code{\link{BDMCMC}} and \code{\link{Phat}} 
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
# Posterior probability for all possible links in the graph
round(Phat(output),2)
}
