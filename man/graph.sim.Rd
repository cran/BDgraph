\name{graph.sim}
\alias{graph.sim}

\title{ Graph simulation }

\description{
Simulating undirected graph structures, including 
"\code{random}", "\code{cluster}", "\code{scale-free}", "\code{lattice}", "\code{hub}", "\code{star}", and "\code{circle}".
}

\usage{
graph.sim( p = 10, graph = "random", prob = 0.2, size = NULL, class = NULL, vis = FALSE,
           rewire = 0.05 )
}

\arguments{
  \item{p}{number of variables (nodes).}
  \item{graph}{ undirected graph with options 
                "\code{random}", "\code{cluster}", "\code{smallworld}", "\code{scale-free}", "\code{lattice}", "\code{hub}", "\code{star}", and "\code{circle}". 
               It also could be an adjacency matrix corresponding to a graph structure (an upper triangular matrix in which 
               \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). 
    }
  \item{prob}{ if \code{graph} = "\code{random}", it is the probability that a pair of nodes has a link.}
  \item{size}{number of links in the true graph (graph size).}
  \item{class}{ if \code{graph} = "\code{cluster}", it is the number of classes. }
  \item{vis}{visualize the true graph structure.}
  \item{rewire}{rewiring probability for smallworld network. Must be between 0 and 1.}
}

\value{ The adjacency matrix corresponding to the simulated graph structure, as an object with \code{S3} class \code{"graph"}.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Alexander Christensen }

\seealso{ \code{\link{bdgraph.sim}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
# Generating a 'hub' graph 
adj <- graph.sim( p = 8, graph = "scale-free" )

plot( adj )

adj
}

\keyword{ simulation }
\keyword{graphs}
