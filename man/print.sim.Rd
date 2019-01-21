\name{print.sim}
\alias{print.sim}

\title{ Print function for \code{S3} class \code{"sim"} }

\description{ Prints the information about the type of data, the sample size, the graph type, the number of nodes, number of links and sparsity of the true graph. }

\usage{ \method{print}{sim}( x, ... ) }

\arguments{
  \item{x}{ An object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}.  }
  \item{\dots}{ System reserved (no specific usage). }
}

\references{
Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and Wit, E. C. (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}, Adrian Dobra, and Ernst Wit }

\seealso{ \code{\link{graph.sim}}, \code{\link{bdgraph.sim}} }

\examples{
# Generating multivariate normal data from a 'random' graph 
data.sim <- bdgraph.sim( n = 20, p = 10, vis = TRUE )
   
print( data.sim )
}
