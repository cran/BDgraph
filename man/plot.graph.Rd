\name{plot.graph}
\alias{plot.graph}

\title{ Plot function for \code{S3} class \code{"graph"} }

\description{ Visualizes structure of the graph. }
\usage{
\method{plot}{graph}( x, main = NULL, layout = layout.circle, ... )
}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"graph"}, from function \code{\link{graph.sim}}.  }
  
  \item{main}{ Graphical parameter (see plot).}
  
  \item{layout}{The vertex placement algorithm which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{\dots}{System reserved (no specific usage).}
}
\references{
Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and Wit, E. C. (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v2} 

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Mohammadi, A. and Dobra, A. (2017). The \code{R} Package \pkg{BDgraph} for Bayesian Structure Learning in Graphical Models, \emph{ISBA Bulletin}, 24(4):11-16
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{\code{\link{graph.sim}}, \code{\link{bdgraph.sim}}}

\examples{
# Generating a 'random' graph 
adj <- graph.sim( p = 10, graph = "random" )
plot( adj )
adj
}
