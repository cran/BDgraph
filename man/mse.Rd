\name{mse}
\alias{mse}

\title{ Graph structure comparison }

\description{ Computes (weighted) mean squared error. }

\usage{ mse(pred, actual, weight = FALSE) }

\arguments{
  \item{pred}{
        adjacency matrix corresponding to an estimated graph. 
        It can be an object with \code{S3} class "\code{bdgraph}" from function \code{\link[BDgraph]{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be an object of \code{S3} class \code{"select"}, from the function \code{\link[huge]{huge.select}} of \code{R} package \code{\link[huge]{huge}}.
        It also can be a \code{list} of above objects for comparing two or more different approaches.
        }

  \item{actual}{
    adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link[BDgraph]{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link[BDgraph]{graph.sim}}.
	}
 
  \item{weight}{ for the case of weighted MSE. }
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
}

\author{Reza Mohammadi \email{a.mohammadi@uva.nl}; Lucas Vogels \email{l.f.o.vogels@uva.nl}}

\seealso{ \code{\link[BDgraph]{compare}}, \code{\link[BDgraph]{auc}}, \code{\link[BDgraph]{bdgraph}}, \code{\link[BDgraph]{bdgraph.mpl}}, \code{\link[BDgraph]{bdgraph.sim}}, \code{\link[BDgraph]{plotroc}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim(n = 50, p = 6, size = 7, vis = TRUE)
    
# Running sampling algorithm based on GGMs 
sample.ggm <- bdgraph(data = data.sim, method = "ggm", iter = 10000)
   
# To compute the value of MSE
mse(pred = sample.ggm, actual = data.sim)

# To compute the value of weighted MSE
mse(pred = sample.ggm, actual = data.sim, weight = 0.5)

}
}

\keyword{structure learning}

