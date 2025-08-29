\name{adj2link}
\alias{adj2link}

\title{ Extract links from an adjacency matrix }

\description{ Extract links from an adjacency matrix or an object of calsses "\code{sim}" from function \code{\link[BDgraph]{bdgraph.sim}} and
 \code{"graph"} from function \code{\link[BDgraph]{graph.sim}}.}

\usage{ adj2link(adj) }

\arguments{
  \item{adj}{     
  adjacency matrix corresponding to a graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link[BDgraph]{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link[BDgraph]{graph.sim}}.
    }
}

\value{ matrix corresponding to the extracted links from graph structure.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link[BDgraph]{link2adj}}, \code{\link[BDgraph]{graph.sim}} }

\examples{
# Generating a 'random' graph 
adj <- graph.sim(p = 6, vis = TRUE)

adj2link(adj)
}

\keyword{graphs}
