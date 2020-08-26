\name{plinks}
\alias{plinks}

\title{ Estimated posterior link probabilities }

\description{ 
Provides the estimated posterior link probabilities for all possible links in the graph.
}

\usage{ plinks( bdgraph.obj, round = 2, burnin = NULL ) }

\arguments{
  \item{bdgraph.obj}{An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
  }

  \item{round}{ A value for rounding all probabilities to the specified number of decimal places. }
  \item{burnin}{ The number of burn-in iteration to scape.  }
}

\value{ An upper triangular matrix which corresponds the estimated posterior probabilities for all possible links.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{\code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'circle' graph
data.sim <- bdgraph.sim( n = 70, p = 6, graph = "circle", vis = TRUE )

bdgraph.obj   <- bdgraph( data = data.sim, iter = 10000 )

plinks( bdgraph.obj, round = 2 )
}
}

\keyword{structure learning}
