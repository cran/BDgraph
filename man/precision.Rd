\name{precision}
\alias{precision}

\title{ Estimated precision matrix }

\description{ 
Provides the estimated precision matrix.
}

\usage{ precision( bdgraph.obj, round = 2 ) }

\arguments{
  \item{bdgraph.obj}{object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
  }

  \item{round}{ value for rounding all probabilities to the specified number of decimal places. }
}

\value{ matrix which corresponds the estimated precision matrix.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{\code{\link{bdgraph}}, \code{\link{covariance}}, \code{\link{plinks}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 70, p = 6, graph = "circle", vis = TRUE )

bdgraph.obj   <- bdgraph( data = data.sim )

precision( bdgraph.obj ) # Estimated precision matrix
  
data.sim $ K   # True precision matrix
}
}

\keyword{parameter learning}
