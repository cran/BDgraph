\name{sparsity}
\alias{sparsity}

\title{ Compute the sparsity of a graph }

\description{ Compute the sparsity of a graph (or a network) or an object of calsses \code{"sim"} from function \code{\link{bdgraph.sim}} and
 \code{"graph"} from function \code{\link{graph.sim}}.}

\usage{ sparsity( adj ) }

\arguments{
  \item{adj}{     
  An adjacency matrix corresponding to a graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class \code{"sim"} from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
    }
}

\value{ A value corresponding to the sparsity of a graph.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{graph.sim}}, \code{\link{adj2link}}, \code{\link{link2adj}} }

\examples{
# Generating a 'random' graph 
adj <- graph.sim( p = 10, vis = TRUE )

sparsity( adj )
}

\keyword{graphs}
