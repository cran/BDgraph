\name{graph.sim}
\alias{graph.sim}

\title{ Graph simulation }

\description{
Simulating undirected graph structures, including 
\code{"random"}, \code{"cluster"}, \code{"scale-free"}, \code{"lattice"}, \code{"hub"}, \code{"star"}, and \code{"circle"}.
}

\usage{
graph.sim( p = 10, graph = "random", prob = 0.2, size = NULL, class = NULL, vis = FALSE )
}

\arguments{
  \item{p}{The number of variables (nodes).}
  \item{graph}{ The undirected graph with options 
                \code{"random"}, \code{"cluster"}, \code{"scale-free"}, \code{"lattice"}, \code{"hub"}, \code{"star"}, and \code{"circle"}. 
               It also could be an adjacency matrix corresponding to a graph structure (an upper triangular matrix in which 
               \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). 
    }
  \item{prob}{ If \code{graph="random"}, it is the probability that a pair of nodes has a link.}
  \item{size}{The number of links in the true graph (graph size).}
  \item{class}{ If \code{graph="cluster"}, it is the number of classes. }
  \item{vis}{Visualize the true graph structure.}
}

\value{ The adjacency matrix corresponding to the simulated graph structure, as an object with \code{S3} class \code{"graph"}.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845

Pensar, J. et al (2017) Marginal pseudo-likelihood learning of discrete Markov network structures, \emph{Bayesian Analysis}, 12(4):1195-215
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{bdgraph.sim}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
# Generating a 'hub' graph 
adj <- graph.sim( p = 8, graph = "scale-free" )

plot( adj )

adj
}

\keyword{ simulation }
\keyword{graphs}
