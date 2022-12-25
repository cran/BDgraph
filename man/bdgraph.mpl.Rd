\name{bdgraph.mpl}
\alias{bdgraph.mpl}
        
\title{	Search algorithm in graphical models using marginal pseudo-likehlihood }
\description{
	This function consists of several sampling algorithms for Bayesian model determination in undirected graphical models based on mariginal pseudo-likelihood. % based on birth-death MCMC method. 
	To speed up the computations, the birth-death MCMC sampling algorithms are implemented in parallel using \pkg{OpenMP} in \code{C++}.
}
\usage{
bdgraph.mpl( data, n = NULL, method = "ggm", transfer = TRUE, 
             algorithm = "bdmcmc", iter = 5000, burnin = iter / 2, 
             g.prior = 0.2, g.start = "empty", 
             jump = NULL, alpha = 0.5, save = FALSE, 
             cores = NULL, operator = "or", verbose = TRUE )
}

\arguments{
	\item{data}{
	there are two options: (1) an (\eqn{n \times p}) \code{matrix} or a \code{data.frame} corresponding to the data, 
	(2) an (\eqn{p \times p}) covariance matrix as \eqn{S=X'X} which \eqn{X} is the data matrix 
	(\eqn{n} is the sample size and \eqn{p} is the number of variables). 
	It also could be an object of class "\code{sim}", from function \code{\link{bdgraph.sim}}.
	The input matrix is automatically identified by checking the symmetry.
	}

	\item{n}{number of observations. It is needed if the "\code{data}" is a covariance matrix.}
	  
	\item{method}{
		character with two options "\code{ggm}" (default), "\code{dgm}" and "\code{dgm-binary}". 
		Option "\code{ggm}" is for Gaussian graphical models based on Gaussianity assumption.
		Option "\code{dgm}" is for discrete graphical models for the count data.
		Option "\code{dgm-binary}" is for discrete graphical models for the data that are binary.
	}

	\item{transfer}{ for only \code{'count'} data which \code{method} = "\code{dgm}" or \code{method} = "\code{dgm-binary}". }

	\item{algorithm}{
		character with two options "\code{bdmcmc}" (default) and "\code{rjmcmc}". 
		Option "\code{bdmcmc}" is based on birth-death MCMC algorithm.
		Option "\code{rjmcmc}" is based on reverible jump MCMC algorithm.
		Option \code{"hc"}     is based on hill-climbing algorithm; this algorithm is only for count data which \code{method} = "\code{dgm}" or \code{method} = "\code{dgm-binary}".
	}
	
	\item{iter}{   number of iteration for the sampling algorithm. }
	\item{burnin}{ number of burn-in iteration for the sampling algorithm. }

	\item{g.prior}{
		for determining the prior distribution of each edge in the graph. 
		There are two options: a single value between \eqn{0} and \eqn{1} (e.g. \eqn{0.5} as a noninformative prior) 
		or an (\eqn{p \times p}) matrix with elements between \eqn{0} and \eqn{1}.
	}
	
	\item{g.start}{
		corresponds to a starting point of the graph. It could be an (\eqn{p \times p}) matrix, "\code{empty}" (default), or "\code{full}". 
		Option "\code{empty}" means the initial graph is an empty graph and "\code{full}" means a full graph. 
		It also could be an object with \code{S3} class "\code{bdgraph}" of \code{R} package \code{\link[BDgraph]{BDgraph}} or the class \code{"ssgraph"} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}; 
		this option can be used to run the sampling algorithm from the last objects of previous run (see examples).     
	} 

	\item{jump}{
		it is only for the BDMCMC algorithm (\code{algorithm} = "\code{bdmcmc}").
		It is for simultaneously updating multiple links at the same time to update graph in the BDMCMC algorithm.
	}	
	\item{alpha}{ value of the hyper parameter of Dirichlet, which is a prior distribution. }
	
	\item{save}{
		logical: if FALSE (default), the adjacency matrices are NOT saved. 
		If TRUE, the adjacency matrices after burn-in are saved.
	}
	\item{cores}{ number of cores to use for parallel execution. 
	   The case \code{cores} = "\code{all}" means all CPU cores to use for parallel execution. 
	   %The default is to use "all" CPU cores of the computer.
	   }	
	\item{operator}{ character with two options "\code{or}" (default) and  "\code{and}". It is for hill-climbing algorithm. }	   	
	\item{verbose}{ logical: if TRUE (default), report/print the MCMC running time. }	
}

\value{
	An object with \code{S3} class "\code{bdgraph}" is returned:
	
	\item{p_links}{ upper triangular matrix which corresponds the estimated posterior probabilities of all possible links. }
	
	For the case "\code{save} = \code{TRUE}" is returned:

	\item{sample_graphs}{ vector of strings which includes the adjacency matrices of visited graphs after burn-in.}
	\item{graph_weights}{ vector which includes the waiting times of visited graphs after burn-in. }

	\item{all_graphs}{vector which includes the identity of the adjacency matrices for all iterations after burn-in. 
		  It is needed for monitoring the convergence of the BD-MCMC algorithm.}

	\item{all_weights}{vector which includes the waiting times for all iterations after burn-in. It is needed for monitoring the convergence of the BD-MCMC algorithm.}
}

\references{
Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215, \doi{10.1214/16-BA1032}

Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}, Adrian Dobra, and Johan Pensar }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.dw}}, \code{\link{bdgraph.sim}}, \code{\link{summary.bdgraph}}, \code{\link{compare}} }

\examples{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 70, p = 5, size = 7, vis = TRUE )
   
bdgraph.obj <- bdgraph.mpl( data = data.sim, iter = 500 )
  
summary( bdgraph.obj )
   
# To compare the result with true graph
compare( bdgraph.obj, data.sim, main = c( "Target", "BDgraph" ) )
}

\keyword{sampling algorithms}
\keyword{structure learning}
\keyword{iteration}
