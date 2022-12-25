\name{bdgraph.dw}
\alias{bdgraph.dw}

\title{	Search algorithm for Gaussian copula graphical models for count data }
\description{
	This function consists of several sampling algorithms for Bayesian structure learning in undirected graphical models for count data. It is based on Gaussian copula graphical models with discrete Weibull distributed marginals.
	To speed up the computations, the birth-death MCMC sampling algorithms are implemented in parallel using \pkg{OpenMP} in \code{C++}.
}
\usage{
bdgraph.dw( data, x = NULL, formula = y ~ ., 
            n = NULL, algorithm = "bdmcmc", iter = 5000, 
            burnin = iter / 2, g.prior = 0.2, df.prior = 3,
            ZI = FALSE, iter_bdw = 5000,
            g.start = "empty", jump = NULL, save = FALSE, 
            q = NULL, beta = NULL, pii = NULL,
            cores = NULL, threshold = 1e-8, verbose = TRUE )
}

\arguments{
    \item{data}{
	    (\eqn{n \times p}) \code{matrix} or a \code{data.frame} corresponding to the data on the \code{p} nodes of the graph. It can also be an object of class "\code{sim}", from the function \code{\link{bdgraph.sim}}.
	    }

	\item{x}{ (\eqn{n \times k}) \code{matrix} or a \code{data.frame} corresponding to the predictors. }

	\item{formula}{ object of class \link{formula} as a symbolic description of the model for linking each node to the predictors. For the case of \code{data.frame}, it is taken as the model frame (see \code{\link{model.frame})}.}

\item{n}{number of observations. It is needed if the "\code{data}" is a covariance matrix.}
	  
	\item{algorithm}{
		character with two options "\code{bdmcmc}" (default) and "\code{rjmcmc}". 
		Option "\code{bdmcmc}" is based on a birth-death MCMC algorithm.
		Option "\code{rjmcmc}" is based on a reversible jump MCMC algorithm.
	}
	
	\item{iter}{number of iterations for the sampling algorithm for graph learning.}
	
	\item{burnin}{number of burn-in iterations for the sampling algorithm for graph learning.}

	\item{g.prior}{
		for determining the prior distribution of each edge in the graph. 
		There are two options: a single value between \eqn{0} and \eqn{1} (e.g. \eqn{0.5} as a noninformative prior) 
		or a (\eqn{p \times p}) matrix with elements between \eqn{0} and \eqn{1}.
	}
	
	\item{df.prior}{degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}, which is a prior distribution for the precision matrix.}

  \item{ZI}{logical. If FALSE (default), the conditional distribution of each response variable is assumed to be Discrete Weibull given the predictors \code{x}. If TRUE, a zero-inflated model will be applied to each response. \code{ZI} can be passed also as a vector, in order to specify which of the (\code{p} variables) should be fitted with zero-inflation (TRUE) or not (FALSE).
	}
	\item{iter_bdw}{number of iterations for the sampling algorithm to estimate the regression parameters for the Discrete Weibull distribution. It is passed to                     the \code{\link{bdw.reg}} function. }

	\item{g.start}{
		corresponds to a starting point of the graph. It could be an (\eqn{p \times p}) matrix, "\code{empty}" (default), or "\code{full}". 
		Option "\code{empty}" means that the initial graph is an empty graph and "\code{full}" means a full graph. 
		It also could be an object with \code{S3} class "\code{bdgraph}" of \code{R} package \code{\link[BDgraph]{BDgraph}} or the class \code{"ssgraph"} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}; 
		this option can be used to run the sampling algorithm from the last objects of the previous run (see examples).     
	} 

	\item{jump}{
		it is only for the BDMCMC algorithm (\code{algorithm} = "\code{bdmcmc}").
		It is for simultaneously updating multiple links at the same time while updating the graph in the BDMCMC algorithm.
	}	
	\item{save}{
		logical: if FALSE (default), the adjacency matrices are NOT saved. 
		If TRUE, the adjacency matrices  after burn-in are saved.
	}
	
    \item{q, beta}{ parameters of the discrete Weibull distribution used for the marginals. They should be given either as a (\eqn{n \times p}) \code{matrix} (if covariates are present) or as a vector (if covariates are not present). If NULL (default), these parameters are estimated by the \code{\link{bdw.reg}} function.}
          
    \item{pii}{ vector of zero-inflation parameters of the zero-inflated discrete Weibull distributions used for the marginals. If NULL (default), this parameter is estimated by the \code{\link{bdw.reg}} function when \code{ZI = TRUE}. }

	\item{cores}{ number of cores to use for parallel execution. 
	   The case \code{cores} = "\code{all}" means all CPU cores to use for parallel execution. 
	   %The default is to use "all" CPU cores of the computer.
    }	
    \item{threshold}{ threshold value for the convergence of the sampling algorithm from G-Wishart for the precision matrix.}
	\item{verbose}{ logical: if TRUE (default), report/print the MCMC running time. }	
}

\value{
	An object with \code{S3} class "\code{bdgraph}" is returned, containing:
	
	\item{p_links}{ upper triangular matrix corresponding to the estimated posterior probabilities of all possible links. }
	
	\item{K_hat}{ posterior estimation of the precision matrix. }
	
	\item{sample_marginals}{ posterior samples of the regression coefficients of the marginal distributions. }

	For the case "\code{save} = \code{TRUE}", the code returns:

	\item{sample_graphs}{ vector of strings which includes the adjacency matrices of the graphs visited after burn-in.}
	\item{graph_weights}{ vector which includes the waiting times of the graphs visited after burn-in. }

	\item{all_graphs}{vector which includes the identity of the adjacency matrices for all iterations after burn-in. 
		  It is needed for monitoring the convergence of the BDMCMC algorithm.}

	\item{all_weights}{vector which includes the waiting times for all iterations after burn-in. 
	                   It is needed for monitoring the convergence of the BDMCMC algorithm.}
}

\references{
Vinciotti, V., Behrouzi, P., and Mohammadi, R. (2022) Bayesian structural learning of microbiota systems from count metagenomic data, \emph{arXiv preprint}, \doi{10.48550/arXiv.2203.10118}

Peluso, A., Vinciotti, V., and Yu, K. (2018) Discrete Weibull generalized additive model: an application to count fertility, \emph{Journal of the Royal Statistical Society: Series C}, 68(3):565-583, \doi{10.1111/rssc.12311}  

Haselimashhadi, H., Vinciotti, V., and Yu, K. (2018) A novel Bayesian regression model for counts with an application to health data, \emph{Journal of Applied Statistics,} 45(6):1085-1105, \doi{10.1080/02664763.2017.1342782}

Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03}  

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}, Veronica Vinciotti, and Pariya Behrouzi }

\seealso{ \code{\link{bdgraph}}, \link{bdgraph.mpl}, \link{bdw.reg}, \code{\link{bdgraph.sim}}, \code{\link{summary.bdgraph}}, \code{\link{compare}} }

\examples{
\dontrun{
# - - Example 1

# Generating multivariate Discrete Weibull data based on 'random' graph
data.sim <- bdgraph.sim( n = 100, p = 10, type = "dw", vis = TRUE )

bdgraph.obj <- bdgraph.dw( data = data.sim, iter = 5000 )

summary( bdgraph.obj )

# To compare the result with true graph
compare( bdgraph.obj, data.sim, main = c( "Target", "BDgraph" ), vis = TRUE )

# - - Example 2

# Generating multivariate Discrete Weibull data based on a 'scale-free' graph
data.sim <- bdgraph.sim( n = 100, p = 10, type = "dw", graph = "scale-free", vis = TRUE )

bdgraph.obj <- bdgraph.dw( data = data.sim, iter = 10000 )

summary( bdgraph.obj )

compare( bdgraph.obj, data.sim, main = c( "Target", "BDgraph" ), vis = TRUE )
}	  
}

\keyword{sampling algorithms}
\keyword{structure learning}
\keyword{iteration}
\keyword{Discrete Weibull}
