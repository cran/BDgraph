\name{pgraph}
\alias{pgraph}

\title{ Posterior probabilities of the graphs }

\description{
Provides the estimated posterior probabilities for the most likely graphs or a specific graph.
}

\usage{ pgraph( bdgraph.obj, number.g = 4, adj = NULL ) }

\arguments{
  \item{bdgraph.obj}{ object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. }
  \item{number.g}{ number of graphs with the highest posterior probabilities to be shown.
            This option is ignored if \code{'adj'} is specified. }
  \item{adj}{ adjacency matrix corresponding to a graph structure. It is an upper triangular matrix in which 
            \eqn{a_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{a_{ij}=0}.
            It also can be an object of \code{S3} class "\code{sim}", from function \code{\link{bdgraph.sim}}. }
}

\value{
\item{ selected_g }{ adjacency matrices which corresponding to the graphs with the highest posterior probabilities.}

\item{ prob_g }{vector of the posterior probabilities of the graphs corresponding to \code{'selected\_g'}.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 6, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim, save = TRUE )
   
# Estimated posterior probability of the true graph
pgraph( bdgraph.obj, adj = data.sim )
   
# Estimated posterior probability of first and second graphs with highest probabilities
pgraph( bdgraph.obj, number.g = 2 )
}
}

\keyword{structure learning}
