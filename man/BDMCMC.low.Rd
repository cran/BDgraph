\name{bdmcmc.low}
\alias{bdmcmc.low}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BDMCMC algorithm for low-dimensional graphs
}
\description{
This function is the BDMCMC algorithm for model selection in undirected Gaussian graphical
models based on birth-death MCMC methodology. The algorithm is not fast, so
it is suitable only for low-dimensional graphs (roughly graphical models with less 
than 8 nodes).
}
\usage{
bdmcmc.low(data, n = NULL, meanzero = FALSE, iter = 5000, burn = floor(iter/2),
           skip = 1, gamma.b = 1, prior.g = "Uniform", b = 3, D = NULL, 
		           A = "full", MCiter = 10, summary = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
it could be a matrix or a data.frame of the data ( \eqn{p \times n}{p x n} matrix) or it could be a
covariance matrix as \eqn{S=x'x}{S = x' x} which \eqn{x}{x} is a matrix of the data.
}
  \item{n}{
the number of observations.
}
  \item{meanzero}{
logical: if TRUE that means the data have zero mean; if FALSE (default) the mean of the data should not be zero.
}
  \item{iter}{
the number of iterations for the BDMCMC algorithm.
}
  \item{burn}{
the number of burn-in iterations for the BDMCMC algorithm.
}
  \item{skip}{
it is the option for regularly saving part of the iterations. Default is 1 (for saving all iterations).
}
  \item{gamma.b}{
the birth rates for the birth-death process.
}
  \item{prior.g}{
a character for selecting a prior distribution for the graph. It can be either "Uniform" (default) or "Poisson".
"Uniform" means discrete uniform distribution for prior distribution of the graph 
(\eqn{G \sim DU(\mathcal{G})}{G ~ DU(G)} in which \eqn{\mathcal{G}}{G} is all possible graphical models).
Also, "Poisson" means prior distribution for number of edges in the graph G has a truncated Poisson distribution
(\eqn{degree(G) \sim TP(\gamma_b)}{degree(G) ~ TP(gamma.b)} in which for our model \eqn{\gamma_b}{gammab} equal 
with a 'gamma.b' value).
}
  \item{b}{
the degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}{WG(b,D)}. G-Wishart distribution is the prior
distribution of precision matrix.
}
  \item{D}{
the positive definite matrix for G-Wishart distribution, \eqn{W_G(b,D)}{WG(b,D)}. G-Wishart distribution is the prior 
distribution of precision matrix.
}
  \item{A}{
matrix A shows the starting graph for BDMCMC algorithm. It is an upper triangular matrix in which 
\eqn{a_{ij}=1}{aij = 1} if there is a link between nodes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}.
It can be either "full" (default) or "emply" or "glasso". "full" means we start with full graph. "emply" means we start
with empty graph. "glasso" means we start with a graph which is suggested by glasso method according to huge package. 
We can use this option for high-dimentional cases which we need less iterations to converge, because we start with the 
graph which is almost near to the true graph.  
} 
  \item{MCiter}{
the number of iterations for Monte Carlo approximation of normalizing constant in G-Wishart distribution. 
We need it only when number of nodes in the graph is less than 8.
}
  \item{summary}{
logical: if TRUE you will see the summary result of BDMCMC algorithm; if FALSE (default), the output will be a 
list as below.
}
  \item{verbose}{
logical: if verbose = TRUE, tracing number of iterations in BDMCMC algorithm printing is abled. The default value is TRUE.
}
}
\value{
Returns an output of the BDMCMC algorithm which is list-like and contains the following:
\item{As}{a list which includes the adjacency matrix for all iterations after burn-in iterations.}
\item{Ks}{a list which includes the precision matrix for all iterations after burn-in iterations.}
\item{lambda}{a vector which includes the waiting times for all iterations after burn-in iterations.}
\item{allA}{a list which includes all adjacency matrix for all iterations (includes burn-in iterations). 
We need list 'allA' for checking the convergency of the BDMCMC algorithm.}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdmcmc}} and \code{\link{bdmcmc.high}}
}
\examples{
\dontrun{
  require(MASS)
  
  p <- 5 # number of nodes
  # "truK" is the precision matrix of true graph
  truK <- diag(p)
  for (i in 1:(p-1)) truK[i,i+1] <- truK[i+1,i] <- 0.5
  truK[1,p] <- truK[p,1] <- 0.4
  truK # precision matrix of the true graph
  
  # generate the data (150 observations) from multivariate normal 
  # distribution with mean zero and percision matrix "truK"
  data <- mvrnorm(150, c(rep(0,p)), solve(truK))
  
  bdmcmc.low(data, meanzero = T, iter = 2000, summary = T)
  }
}

