\name{bdgraph}
\alias{bdgraph}

\title{	Search algorithm in graphical models }
\description{
	As the main function of the \pkg{BDgraph} package, 
	this function consists of several MCMC sampling algorithms for Bayesian model determination in undirected graphical models. % based on birth-death MCMC method. 
	To speed up the computations, the birth-death MCMC sampling algorithms are implemented in parallel using \pkg{OpenMP} in \code{C++}.
}
\usage{
bdgraph( data, n = NULL, method = "ggm", algorithm = "bdmcmc", iter = 5000, 
         burnin = iter / 2, not.cont = NULL, g.prior = 0.5, df.prior = 3, 
         g.start = "empty", jump = NULL, save = FALSE, print = 1000, 
         cores = NULL, threshold = 1e-8 )
}

\arguments{
	\item{data}{
	There are two options: (1) an (\eqn{n \times p}) matrix or a \code{data.frame} corresponding to the data, 
	(2) an (\eqn{p \times p}) covariance matrix as \eqn{S=X'X} which \eqn{X} is the data matrix 
	(\eqn{n} is the sample size and \eqn{p} is the number of variables). 
	It also could be an object of class \code{"sim"}, from function \code{\link{bdgraph.sim}}.
	The input matrix is automatically identified by checking the symmetry.
	}

	\item{n}{The number of observations. It is needed if the \code{"data"} is a covariance matrix.}
	  
	\item{method}{
		A character with two options \code{"ggm"} (default) and \code{"gcgm"}. 
		Option \code{"ggm"} is for Gaussian graphical models based on Gaussianity assumption.
		Option \code{"gcgm"} is for Gaussian copula graphical models for the data that not follow Gaussianity assumption (e.g. continuous non-Gaussian, discrete, or mixed dataset).
	}
	
	\item{algorithm}{
		A character with two options \code{"bdmcmc"} (default) and \code{"rjmcmc"}. 
		Option \code{"bdmcmc"} is based on birth-death MCMC algorithm.
		Option \code{"rjmcmc"} is based on reverible jump MCMC algorithm.
	}
	
	\item{iter}{The number of iteration for the sampling algorithm.}
	\item{burnin}{The number of burn-in iteration for the sampling algorithm.}

	\item{not.cont}{ For the case \code{method = "gcgm"}, a vector with binary values in which \eqn{1} indicates not continuous variables. }

	\item{g.prior}{
		For determining the prior distribution of each edge in the graph. 
		There are two options: a single value between \eqn{0} and \eqn{1} (e.g. \eqn{0.5} as a noninformative prior) 
		or an (\eqn{p \times p}) matrix with elements between \eqn{0} and \eqn{1}.
	}
	
	\item{df.prior}{
		The degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}, which is a prior distribution of the precision matrix.
	}

	\item{g.start}{
		Corresponds to a starting point of the graph. It could be an (\eqn{p \times p}) matrix, \code{"empty"} (default), or \code{"full"}. 
		Option \code{"empty"} means the initial graph is an empty graph and \code{"full"} means a full graph. 
		It also could be an object with \code{S3} class \code{"bdgraph"} of \code{R} package \code{\link[BDgraph]{BDgraph}} or the class \code{"ssgraph"} of \code{R} package \code{\link[ssgraph]{ssgraph}}; 
		this option can be used to run the sampling algorithm from the last objects of previous run (see examples).     
	} 

	\item{jump}{
		It is only for the BDMCMC algorithm (\code{algorithm = "bdmcmc"}).
		It is for simultaneously updating multiple links at the same time to update graph in the BDMCMC algorithm.
	}	
	\item{save}{
		Logical: if FALSE (default), the adjacency matrices are NOT saved. 
		If TRUE, the adjacency matrices after burn-in are saved.
	}
	\item{print}{ Value to see the number of iteration for the MCMC algorithm. }	
	\item{cores}{ The number of cores to use for parallel execution. 
	   The case \code{cores="all"} means all CPU cores to use for parallel execution. 
	   %The default is to use "all" CPU cores of the computer.
    }	
    \item{threshold}{ The threshold value for the convergence of sampling algorithm from G-Wishart for the precision matrix.}
}

\value{
	An object with \code{S3} class \code{"bdgraph"} is returned:
	
	\item{p_links}{ An upper triangular matrix which corresponds the estimated posterior probabilities of all possible links. }
	
	\item{K_hat}{ The posterior estimation of the precision matrix. }
	
	For the case "save = TRUE" is returned:

	\item{sample_graphs}{ A vector of strings which includes the adjacency matrices of visited graphs after burn-in.}
	\item{graph_weights}{ A vector which includes the waiting times of visited graphs after burn-in. }

	\item{all_graphs}{A vector which includes the identity of the adjacency matrices for all iterations after burn-in. 
		  It is needed for monitoring the convergence of the BD-MCMC algorithm.}

	\item{all_weights}{A vector which includes the waiting times for all iterations after burn-in. 
	                   It is needed for monitoring the convergence of the BD-MCMC algorithm.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. and Dobra A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph.mpl}}, \code{\link{bdgraph.sim}}, \code{\link{summary.bdgraph}}, \code{\link{compare}} }

\examples{
\dontrun{
# --- Example 1
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 20, p = 6, size = 7, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, iter = 1000 )
  
summary( bdgraph.obj )
   
# To compare our result with true graph
compare( data.sim, bdgraph.obj, main = c( "Target", "BDgraph" ) )
   
# Running algorithm with starting points from previous run
bdgraph.obj2 <- bdgraph( data = data.sim, g.start = bdgraph.obj )
    
compare( data.sim, bdgraph.obj, bdgraph.obj2, 
         main = c( "Target", "Frist run", "Second run" ) )

# --- Example 2
# Generating mixed data from a 'scale-free' graph
data.sim <- bdgraph.sim( n = 50, p = 6, type = "mixed", graph = "scale-free", vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )
  
summary( bdgraph.obj )
   
compare( data.sim, bdgraph.obj )
}	  
}
