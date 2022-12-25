\name{plot.graph}
\alias{plot.graph}

\title{ Plot function for \code{S3} class \code{"graph"} }

\description{ Visualizes structure of the graph. }
\usage{
\method{plot}{graph}( x, cut = 0.5, mode = "undirected", diag = FALSE, main = NULL, 
           layout = igraph::layout_with_fr, vertex.size = 2, vertex.color = "orange", 
           vertex.frame.color = "orange", vertex.label = NULL, vertex.label.dist = 0.5, 
           vertex.label.color = "blue", edge.color = "lightblue", ... )
}

\arguments{
  \item{x}{   object of \code{S3} class \code{"graph"}, from function \code{\link{graph.sim}}.  }
  \item{cut}{ for the case where input '\code{x}' is the object of class "bdgraph" or "ssgraph". Threshold for including the links in the selected graph based on the estimated posterior probabilities of the links. } 
  
  \item{mode}{ type of graph which is according to \code{R} package \code{\link[igraph]{igraph}}.}
  \item{diag}{ logical which is according to \code{R} package \code{\link[igraph]{igraph}}.}
  \item{main}{ graphical parameter (see plot).}

  \item{layout}{ vertex placement which is according to \code{R} package \code{\link[igraph]{igraph}}; For different layouts, see \code{\link[igraph]{layout}} of \code{R} package \code{igraph}. }
  \item{vertex.size}{ vertex size which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  \item{vertex.color}{ vertex color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  \item{vertex.frame.color}{ vertex frame color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{vertex.label}{ vertex label. The default vertex labels are the vertex ids. }
  \item{vertex.label.dist}{ vertex label distance which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  \item{vertex.label.color}{ vertex label color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{edge.color}{ edge color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{\dots}{ additional plotting parameters. For the complete list, see \code{\link[igraph]{igraph.plotting}} of \code{R} package \code{igraph}. }
}
\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{graph.sim}}, \code{\link{bdgraph.sim}}, \code{\link{plot.igraph}} }

\examples{
# Generating a 'scale-free' graph 
adj <- graph.sim( p = 20, graph = "scale-free" )

plot( adj )
}

\keyword{hplot}

