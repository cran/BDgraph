\name{bdgraph}
\alias{bdgraph}

\title{
Graph selection based on birth-death MCMC algorithm 
}
\description{
This function is the BD-MCMC algorithm for model selection in undirected Gaussian 
graphical models based on birth-death MCMC methodology. 
}
\usage{
bdgraph( data, n = NULL, method = "exact", npn = "normal", iter = 5000, 
		 burnin = floor(iter / 2), b = 3, D = NULL, Gstart = "empty", 
		 Kstart = NULL, trace = TRUE )
}

\arguments{

\item{data}{
It could be a matrix or a data.frame of the data ( \eqn{n \times p}) or it could be a
covariance matrix as \eqn{S=x'x} which \eqn{x} is the data matrix. It also could be an object of class "simulate".
}

\item{n}{The number of observations. It is needed if the "data" is a covariance matrix.}
  
\item{method}{
A character with three options "exact" (defult), "approx" and "copula" to determine a type of BD-MCMC algorithm. 
Option "exact" is a exact sampling from joint posterior distributions, without any approximation.
Option "approx" is an approximation approach for sampling from joint posterior distributions. 
See first reference for more details.
Option "copula" is for non-Gaussian, discrete, or mixed data sets.
}

\item{npn}{ 
A character with four options "normal" (default), "shrinkage", "truncation", and "skeptic". 
Option "normal" means data are coming from multivariate normal distribution. 
Option "shrinkage" is for the shrunken transformation, option "truncation" is for 
the truncated transformation, option "skeptic" is for the non-paranormal skeptic transformation. 
See \code{\link{bdgraph.npn}}.
}

\item{iter}{The number of iteration for the BD-MCMC algorithm.}
\item{burnin}{The number of burn-in iteration for the BD-MCMC algorithm.}

\item{b}{
The degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}. G-Wishart distribution is the prior
distribution of precision matrix. The default is 3.
}
\item{D}{
The positive definite matrix for G-Wishart distribution, \eqn{W_G(b,D)}. G-Wishart distribution is the prior 
distribution of precision matrix. The default is an identity matrix.
}
\item{Gstart}{
It is corresponding to a starting point for graph. It can be either "full" (default), "empty", or an object with S3 class "bdgraph". 
Option "full" means the initial graph is a full graph and "empty" means a empty graph. 
Gstart also could be an object with S3 class "bdgraph"; With this option we could run the bdgraph algorithm from last objects of 
previous bdgraph run (see examples).     
} 

\item{Kstart}{It is corresponding to a starting point for precision matrix. It has to be a positive definite matrix.}

\item{trace}{Logical: if TRUE (default), prints out algorithm progress. }

}

\value{
an object with S3 class "bdgraph" as follow:

\item{sampleGraphs}{A vector which includes the adjacency matrices for all iteration after burn-in.}
\item{graphWeights}{A vector which includes the waiting times for all iteration after burn-in.}

\item{allGraphs}{A vector which includes the adjacency matrices for all iteration (includes burn-in iteration). 
      It is needed to check the convergency of the BD-MCMC algorithm.}

\item{allWeights}{A vector which includes the waiting times for all iteration (includes burn-in iteration).}

\item{Khat}{Estimation for precision matrix which is a mean of all samples from precision matrices.}
}

\references{
Mohammadi, A. and E. C. Wit (2014). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{ \code{\link{bdgraph.sim}} }

\examples{
\dontrun{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 20, p = 6, size = 7, vis = TRUE )
  
  output <- bdgraph( data = data.sim, iter = 1000 )
  
  summary(output)
  
  # To compare our result with true graph
  compare( data.sim, output, colnames = c("True graph", "BDgraph") )
  
  output2 <- bdgraph( data = data.sim, iter = 5000, g.start = output )
  
  compare( data.sim, output, output2, colnames = c("True graph", "Frist run", "Second run") )
  }
}

