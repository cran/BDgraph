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
         burnin = iter / 2, not.cont = NULL, g.prior = 0.2, df.prior = 3, 
         g.start = "empty", jump = NULL, save = FALSE, 
         cores = NULL, threshold = 1e-8, verbose = TRUE, nu = 1 )
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
		character with two options "\code{ggm}" (default) and "\code{gcgm}". 
		Option "\code{ggm}" is for Gaussian graphical models based on Gaussianity assumption.
		Option "\code{gcgm}" is for Gaussian copula graphical models for the data that not follow Gaussianity assumption (e.g. continuous non-Gaussian, count, or mixed dataset).
	}
	
	\item{algorithm}{
		character with two options "\code{bdmcmc}" (default) and "\code{rjmcmc}". 
		Option "\code{bdmcmc}" is based on birth-death MCMC algorithm.
		Option "\code{rjmcmc}" is based on reverible jump MCMC algorithm.
		Option "\code{bd-dmh}" is based on birth-death MCMC algorithm using double Metropolis Hasting.
		Option "\code{rj-dmh}" is based on reverible jump MCMC algorithm using double Metropolis Hasting.
	}
	
	\item{iter}{number of iteration for the sampling algorithm.}
	\item{burnin}{number of burn-in iteration for the sampling algorithm.}

	\item{not.cont}{ for the case \code{method = "gcgm"}, a vector with binary values in which \eqn{1} indicates not continuous variables. }

	\item{g.prior}{
		for determining the prior distribution of each edge in the graph. 
		There are two options: a single value between \eqn{0} and \eqn{1} (e.g. \eqn{0.5} as a noninformative prior) 
		or an (\eqn{p \times p}) matrix with elements between \eqn{0} and \eqn{1}.
	}
	
	\item{df.prior}{
		degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}, which is a prior distribution of the precision matrix.
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
	\item{save}{
		logical: if FALSE (default), the adjacency matrices are NOT saved. 
		If TRUE, the adjacency matrices after burn-in are saved.
	}
	\item{cores}{ number of cores to use for parallel execution. 
	   The case \code{cores} = "\code{all}" means all CPU cores to use for parallel execution. 
	   %The default is to use "all" CPU cores of the computer.
    }	
    \item{threshold}{ threshold value for the convergence of sampling algorithm from G-Wishart for the precision matrix.}
	\item{verbose}{ logical: if TRUE (default), report/print the MCMC running time. }
	\item{nu}{ prior parameter for option \code{method} = "\code{tgm}".}
}

\value{
	An object with \code{S3} class "\code{bdgraph}" is returned:
	
	\item{p_links}{ upper triangular matrix which corresponds the estimated posterior probabilities of all possible links. }
	
	\item{K_hat}{ posterior estimation of the precision matrix. }
	
	For the case "\code{save} = \code{TRUE}" is returned:

	\item{sample_graphs}{ vector of strings which includes the adjacency matrices of visited graphs after burn-in.}
	\item{graph_weights}{ vector which includes the waiting times of visited graphs after burn-in. }

	\item{all_graphs}{vector which includes the identity of the adjacency matrices for all iterations after burn-in. 
		  It is needed for monitoring the convergence of the BD-MCMC algorithm.}

	\item{all_weights}{vector which includes the waiting times for all iterations after burn-in. 
	                   It is needed for monitoring the convergence of the BD-MCMC algorithm.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171}
  
Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. and Dobra A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph.mpl}}, \code{\link{bdgraph.dw}}, \code{\link{bdgraph.sim}}, \code{\link{summary.bdgraph}}, \code{\link{compare}} }

\examples{
\dontrun{
set.seed( 10 )

# - - Example 1

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 100, p = 10, size = 15, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, iter = 1000, save = TRUE )
  
summary( bdgraph.obj )

# Confusion Matrix
conf.mat( actual = data.sim, pred = bdgraph.obj )

conf.mat.plot( actual = data.sim, pred = bdgraph.obj )
   
# To compare our result with true graph
compare( bdgraph.obj, data.sim, main = c( "Target", "BDgraph" ), vis = T )

# Running algorithm with starting points from previous run
bdgraph.obj2 <- bdgraph( data = data.sim, g.start = bdgraph.obj )
    
compare( list( bdgraph.obj, bdgraph.obj2 ), data.sim, 
         main = c( "Target", "Frist run", "Second run" ) )

# - - Example 2

# Generating mixed data from a 'scale-free' graph
data.sim <- bdgraph.sim( n = 200, p = 7, type = "mixed", graph = "scale-free", vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, method = "gcgm" )
  
summary( bdgraph.obj )
   
compare( bdgraph.obj, data.sim, vis = T )

conf.mat( actual = data.sim, pred = bdgraph.obj )

conf.mat.plot( actual = data.sim, pred = bdgraph.obj )
}	  
}

\keyword{sampling algorithms}
\keyword{structure learning}
\keyword{iteration}
