\name{bf}
\alias{bf}

\title{ Bayes factor for two graphs }

\description{ Compute the Bayes factor for two graph structures. }

\usage{ 
    bf( num, den, bdgraph.obj, log = TRUE ) 
}

\arguments{
    \item{num,
        den}{
        adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
        It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
        It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
    }

  \item{bdgraph.obj}{object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
  }
    
    \item{log}{ character value. If TRUE the Bayes factor is given as log(BF). }
}

\value{ single numeric value, the Bayes factor of the two graph structures \code{num} and \code{den}.}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{compare}}, \code{\link{bdgraph.sim}}}

\examples{
    \dontrun{
        # Generating multivariate normal data from a 'circle' graph
        data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )

        # Running sampling algorithm
        bdgraph.obj <- bdgraph( data = data.sim )

        graph_1 <- graph.sim( p = 6, vis = TRUE )
        
        graph_2 <- graph.sim( p = 6, vis = TRUE )

        bf( num = graph_1, den = graph_2, bdgraph.obj = bdgraph.obj )
    }
}

\keyword{structure learning}

