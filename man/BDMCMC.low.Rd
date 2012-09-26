\name{BDMCMC.low}
\alias{BDMCMC.low}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BDMCMC algorithm for low-dimensional graphs
}
\description{
This function is the BDMCMC algorithm for model selection in Gussian grohical
models based on birth-death MCMC methodology. The algorithm is quitly slow, so
it is suitable only for low-dimensional problems (graphical models with less 
than 8 nodes).
}
\usage{
BDMCMC.low(S, n = NULL, iter = 5000, burn = 2000, distance = 1, lambda_b = 1, b = 3, D = NULL, A = NULL, MC.iter = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
It could be our data (matrix p*n) or it could be matrix S=x'x that x is our data.
}
  \item{n}{
Number of observations 
}
  \item{iter}{
Number of iterations for BDMCMC algorithm, default is 5000
}
  \item{burn}{
Number of iterations as burn-in, default is 2000
}
  \item{distance}{
"distance" is the option for regularly saving part of the iterations. Default is 1 (for saving all the iterations).
}
  \item{lambda_b}{
Birth rates for brith-death process, default is 1
}
  \item{b}{
Value for prior distribution of precision matrix, default is 3
}
  \item{D}{
Positive definite matrix for prior distribution of precision matrix, default is identity matrix
}
  \item{A}{
Upper truculer matrix which shows the starting graphs in which a_ij=1 if there 
is a link between notes i and j, otherwise a_ij=0  
}
  \item{MC.iter}{
Number of iterations for Monte Carlo approximation of normalizing constant, default is 10
}
}
\value{
A list with components "Asample", "Ksample", "lambda.s" and "allA", each one with components:
\item{Asample}{A list which includes the adjacency matrix for all iterations after burn_in iterations}
\item{Ksample}{A list which includes the precision matrix for all iterations after burn_in iterations}
\item{lambda_s}{A vector which includes the waiting times for all iterations after burn_in iterations}
\item{allA}{A list which includes all adjacency matrix for all iterations include burn_in iterations to check the conversgency of BDMCMC algorithm}
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
\code{\link{BDMCMC}}, \code{\link{BDMCMC.high}} and \code{\link{Phat}}.
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
\dontrun{output = BDMCMC.low(data) }
# Posterior probability for all possible links in the graph
round(Phat(output),2)

}

