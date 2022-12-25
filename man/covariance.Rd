\name{covariance}
\alias{covariance}

\title{ Estimated covariance matrix }

\description{ 
Provides the estimated covariance matrix.
}

\usage{ covariance( bdgraph.obj, round = 2 ) }

\arguments{
  \item{bdgraph.obj}{object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
  }

  \item{round}{ value for rounding all probabilities to the specified number of decimal places. }
}

\value{ matrix which corresponds the estimated covariance matrix.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{\code{\link{bdgraph}}, \code{\link{precision}}, \code{\link{plinks}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 70, p = 6, graph = "circle", vis = TRUE )

bdgraph.obj <- bdgraph( data = data.sim )

covariance( bdgraph.obj ) # Estimated covariance matrix
  
data.sim $ sigma   # True covariance matrix
}
}

\keyword{parameter learning}
