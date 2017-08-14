\name{plinks}
\alias{plinks}

\title{ Estimated posterior link probabilities }

\description{ 
Provides the estimated posterior link probabilities for all possible links in the graph.
}

\usage{ plinks( bdgraph.obj, round = 2 ) }

\arguments{
  \item{bdgraph.obj}{An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }

  \item{round}{ A value for rounding all probabilities to the specified number of decimal places. }
}

\value{
\item{p_links}{ An upper triangular matrix which corresponds the estimated posterior probabilities for all possible links.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Mohammadi, A., Massam H., and G. Letac (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{\code{\link{bdgraph}} and \code{\link{bdgraph.sim}}}

\examples{
	\dontrun{
	# Generating multivariate normal data from a 'cycle' graph
	data.sim <- bdgraph.sim( n = 70, p = 6, graph = "cycle", vis = TRUE )
   
	bdgraph.obj   <- bdgraph( data = data.sim, iter = 10000 )
  
	plinks( bdgraph.obj, round = 2 )
	}
}

