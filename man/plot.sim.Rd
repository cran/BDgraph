\name{plot.sim}
\alias{plot.sim}

\title{ Plot function for \code{S3} class "\code{sim}" }

\description{ Visualizes structure of the simulated graph for an object of \code{S3} class "\code{sim}", from function \code{\link{bdgraph.sim}}. }
\usage{
\method{plot}{sim}( x, ... )
}

\arguments{
  \item{x}{ object of \code{S3} class "\code{sim}", from function \code{\link{bdgraph.sim}}.  }
  
  \item{\dots}{ additional plotting parameters. See \code{\link[BDgraph]{plot.graph}} and for the complete list \code{\link[igraph]{igraph.plotting}} of \code{R} package \code{igraph}. }
}
\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{graph.sim}}, \code{\link{bdgraph.sim}}, \code{\link{plot.graph}}, \code{\link{plot.igraph}} }

\examples{
set.seed( 10 )

# Generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 10, p = 15 )
  
plot( data.sim )
}

\keyword{hplot}
