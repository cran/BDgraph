\name{link2adj}
\alias{link2adj}

\title{ Extract links from an adjacency matrix }

\description{ Extract links from an adjacency matrix or an object of calsses "\code{sim}" from function \code{\link{bdgraph.sim}} and
 \code{"graph"} from function \code{\link{graph.sim}}.}

\usage{ link2adj( link, p = NULL ) }

\arguments{
  \item{link}{ (\eqn{2 \times p}) \code{matrix} or a \code{data.frame} corresponding to the links from the graph structure. }
  \item{p}{number of nodes of the graph.}
}

\value{ An adjacency matrix corresponding to a graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{adj2link}}, \code{\link{graph.sim}} }

\examples{
# Generating a 'random' graph 
adj <- graph.sim( p = 6, vis = TRUE )

link <- adj2link( adj )

link2adj( link, p = 6 )
   
}

\keyword{graphs}
