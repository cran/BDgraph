\name{bdmcmc}
\alias{bdmcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BDMCMC algorithm for Gaussain graphical model determination
}
\description{
This function is the BDMCMC algorithm for model selection in Gussian grohical
models based on birth-death MCMC methodology. It works for both low 
and high dimentioanl graphs.
}
\usage{
bdmcmc(data, n = NULL, meanzero = FALSE, iter = 5000, burn = 2000, distance = 1, 
gamma.b = 1, Gprior = "Uniform", b = 3, D = NULL, A = NULL, MCiter = 10, print = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
It could be our data (matrix p*n) or could be matrix S=x'x that matrix x is our data.
}
  \item{n}{
The number of observations.
}
  \item{meanzero}{
Logical; if TRUE it is assumed that the data have zero mean; if FALSE (default), the mean of the data should not be zero.
}
  \item{iter}{
The number of iterations for the BDMCMC algorithm.
}
  \item{burn}{
The number of burn-in iterations for the BDMCMC algorithm.
}
  \item{distance}{
It is the option for regularly saving part of the iterations. Default is 1 (for saving all the iterations).
}
  \item{gamma.b}{
Birth rates for the brith-death process.
}
  \item{Gprior}{
A character for selecting a prior distribution for the graph. It can be either "Uniform" (default) or "Poisson".
"Uniform" means discrete uniform distribution for prior distribution of the graph (G~DU(all possible graphs)) 
and "Poisson" means prior distribution for number of links in the graph is truncated Poisson distribution with
rate gamma.b (degree(G)~TP(gamma.b)).
}
  \item{b}{
The degree of freedom for G-Wishart distribution, W_G(b,D). G-Wishart distrubution is the prior distribution 
of precision matrix.
}
  \item{D}{
Positive definite matrix for G-Wishart distribution, W_G(b,D). G-Wishart distrubution is the prior distribution 
of precision matrix.
}
  \item{A}{
Upper truculer matrix in which a_ij=1 if there is a link between notes i and j, otherwise a_ij=0.
This matrix shows the starting graphal model for BDMCMC algorithm.
} 
  \item{MCiter}{
The number of iterations for Monte Carlo approximation of normalizing constant in G-Wishart distribution. 
We need it only when number of nodes in the graph is less than 8.
}
  \item{print}{
logical; if TRUE you will see the print of iteration, time, and number of links in a current state of BDMCMC algorithm in 
each iteration; if FALSE (default), there is no any print order.
}
}
\value{
A list with components "Asample", "Ksample", "lambda.s" and "allA", each one with components:
\item{Asample}{A list which includes the adjacency matrix for all iterations after burn-in iterations.}
\item{Ksample}{A list which includes the precision matrix for all iterations after burn-in iterations.}
\item{lambda.s}{A vector which includes the waiting times for all iterations after burn-in iterations.}
\item{allA}{A list which includes all adjacency matrix for all iterations (include burn-in iterations). 
We need allA list for checking the convergency of the BDMCMC algorithm.}
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
\code{\link{bdmcmc.low}}, \code{\link{bdmcmc.high}} and \code{\link{phat}}
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
  # Posterior edge inclusion probabioities for all possible edges
  phat(output,2)
  }
}

