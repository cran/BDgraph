\name{precision}
\alias{precision}

\title{ Estimated precision matrix }

\description{ 
Provides the estimated precision matrix.
}

\usage{ precision( bdgraph.obj, round = 2 ) }

\arguments{
  \item{bdgraph.obj}{An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph]{ssgraph}} of \code{R} package \code{\link[ssgraph]{ssgraph}}.
  }

  \item{round}{ A value for rounding all probabilities to the specified number of decimal places. }
}

\value{ A matrix which corresponds the estimated precision matrix.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 
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
