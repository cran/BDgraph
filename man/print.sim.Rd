\name{print.sim}
\alias{print.sim}

\title{ Print function for \code{S3} class "\code{sim}" }

\description{ Prints the information about the type of data, the sample size, the graph type, the number of nodes, number of links and sparsity of the true graph. }

\usage{ \method{print}{sim}( x, ... ) }

\arguments{
  \item{x}{ object of \code{S3} class "\code{sim}", from function \code{\link{bdgraph.sim}}.  }
  \item{\dots}{ system reserved (no specific usage). }
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{graph.sim}}, \code{\link{bdgraph.sim}} }

\examples{
# Generating multivariate normal data from a 'random' graph 
data.sim <- bdgraph.sim( n = 20, p = 10, vis = TRUE )
   
print( data.sim )
}

\keyword{print}
