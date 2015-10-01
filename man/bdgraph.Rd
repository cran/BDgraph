\name{bdgraph}
\alias{bdgraph}

\title{	Search algorithm in graphical models }
\description{
	As the main function of the \pkg{BDgraph} package, 
	this function consists of several sampling algorithms for Bayesian model determination in undirected graphical models. % based on birth-death MCMC method. 
}
\usage{
bdgraph( data, n = NULL, method = "ggm", algorithm = "bdmcmc", iter = 5000, 
		         burnin = iter / 2, b = 3, Gstart = "empty", save.all = FALSE )
}

\arguments{
	\item{data}{
	It could be a (\eqn{n \times p}) matrix or a \code{data.frame} of data or a covariance matrix as \eqn{S=X'X} which \eqn{X} is the data matrix. It also could be an object of class \code{"sim"}, from function \code{\link{bdgraph.sim}}.
	}

	\item{n}{The number of observations. It is needed if the \code{"data"} is a covariance matrix.}
	  
	\item{method}{
		A character with two options \code{"ggm"} (defult) and \code{"gcgm"}. 
		Option \code{"ggm"} is for Gaussian graphical models based on Gaussianity assumption.
		Option \code{"gcgm"} is for Gaussian copula graphical models for the data that not follow Gaussianity assumption (e.g. continuous non-Gaussian, discrete, or mixed dataset).
	}

	\item{algorithm}{
		A character with two options \code{"bdmcmc"} (defult) and \code{"rjmcmc"}. 
		Option \code{"bdmcmc"} is based on birth-death MCMC algorithm.
		Option \code{"rjmcmc"} is based on reverible jump MCMC algorithm.
	}
	
	\item{iter}{The number of iteration for the sampling algorithm.}
	\item{burnin}{The number of burn-in iteration for the sampling algorithm.}

	\item{b}{
		The degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}, which is a prior distribution of the precision matrix. The default is 3.
	}
	\item{Gstart}{
		Corresponds to a starting point of the graph. It can be \code{"full"} (default), \code{"empty"}, or an object with \code{S3} class \code{"bdgraph"}. 
		Option \code{"full"} means the initial graph is a full graph and \code{"empty"} means an empty graph. 
		Gstart also could be an object with \code{S3} class \code{"bdgraph"}; With this option we could run the sampling algorithm from last objects of 
		previous run (see examples).     
	} 
	\item{save.all}{
		Logical: if FALSE (default), the adjacency matrices are NOT saved. 
		If TRUE, the adjacency matrices after burn-in are saved.
	}
}

\value{
	An object with \code{S3} class \code{"bdgraph"} is returned:
	
	\item{phat}{ Upper triangular matrix which corresponds the estimated posterior probabilities of all possible links. }
	
	\item{Khat}{ Posterior estimation of the precision matrix. }
	
	For the case "save.all = TRUE" is retrurned:

	\item{sampleGraphs}{ A vector of strings which includes the adjacency matrices of visited graphs after burn-in.}
	\item{graphWeights}{ A vector which includes the waiting times of visited graphs after burn-in. }

	\item{allGraphs}{A vector which includes the identity of the adjacency matrices for all iterations after burn-in. 
		  It is needed for monitoring the convergence of the BD-MCMC algorithm.}

	\item{allWeights}{A vector which includes the waiting times for all iterations after burn-in. It is needed for monitoring the convergence of the BD-MCMC algorithm.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{ \code{\link{bdgraph.sim}}, \code{\link{summary.bdgraph}}, and \code{\link{compare}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 20, p = 6, size = 7, vis = TRUE )
   
output <- bdgraph( data = data.sim, iter = 1000 )
  
summary(output)
   
# To compare our result with true graph
compare( data.sim, output, colnames = c("True graph", "BDgraph") )
   
# Running algorithm with starting points from previous run
output2 <- bdgraph( data = data.sim, iter = 5000, Gstart = output )

compare( data.sim, output, output2, colnames = c("True graph", "Frist run", "Second run") )
   
# Generating mixed data from a 'scale-free' graph
data.sim <- bdgraph.sim( n = 50, p = 6, type = "mixed", graph = "scale-free", vis = TRUE )

output <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )

summary( output )
   
compare( data.sim, output )	  
}
}

