\name{plotcoda}
\alias{plotcoda}

\title{ Convergence plot }

\description{
Visualizes the cumulative occupancy fractions of all possible links in the graph.
It can be used for monitoring the convergence of the sampling algorithms, BDMCMC and RJMCMC.
}

\usage{ plotcoda( bdgraph.obj, thin = NULL, control = TRUE, main = NULL, 
                  verbose = TRUE, ... ) 
}

\arguments{
  \item{bdgraph.obj}{ object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
  }
  
  \item{thin}{ option for getting fast result for a cumulative plot according to part of the iteration.}
    
  \item{control}{ logical: if TRUE (default) and the number of nodes is greater than 15, then 100 links randomly is selected for visualization. } 
  
  \item{main}{ graphical parameter (see plot). }
	\item{verbose}{ logical: if TRUE (default), report/print the calculation progress. }	
  
  \item{\dots}{ system reserved (no specific usage). }  
}

\details{
Note that a spending time for this function depends on the number of nodes. 
%It should be slow for high-dimensional graphs. 
For fast result, you can choose bigger value for the \code{'thin'} option.
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{traceplot}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )
  
bdgraph.obj <- bdgraph( data = data.sim, iter = 10000, burnin = 0 , save = TRUE )
   
plotcoda( bdgraph.obj )  
}
}

\keyword{hplot}
