\name{sparsity}
\alias{sparsity}

\title{ Compute the sparsity of a graph }

\description{ Compute the sparsity of a graph/network or an object of calss \code{"graph"} from function \code{\link{graph.sim}} or an object of calss "\code{sim}" from function \code{\link{bdgraph.sim}}. }

\usage{ sparsity( adj ) }

\arguments{
  \item{adj}{     
  adjacency matrix corresponding to a graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
    }
}

\value{ value corresponding to the graph sparsity which is the proportion of the non-links (non-zero elements) in \code{adj}. }

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{graph.sim}}, \code{\link{adj2link}}, \code{\link{link2adj}} }

\examples{
# Generating a 'random' graph 
adj <- graph.sim( p = 10, graph = "random", prob = 0.4, vis = TRUE )

sparsity( adj )
}

\keyword{graphs}
