\name{print.bdgraph}
\alias{print.bdgraph}

\title{Print function for \code{S3} class \code{"bdgraph"} }

\description{
Prints the information about the selected graph which could be a graph with links for which their estimated posterior probabilities are greater than 0.5
or graph with the highest posterior probability.
It provides adjacency matrix, size and posterior probability of the selected graph.
}

\usage{\method{print}{bdgraph}( x, round = 2, ... )}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{round}{A value to round the probabilities to the specified number of decimal places.}
 
  \item{\dots}{System reserved (no specific usage). }
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Mohammadi, A., Massam H., and G. Letac (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, vis = TRUE )
  
bdgraph.obj <- bdgraph( data = data.sim )
   
print( bdgraph.obj )
}
}
