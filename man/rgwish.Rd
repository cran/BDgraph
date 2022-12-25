\name{rgwish}
\alias{rgwish}

\title{Sampling from G-Wishart distribution}

\description{
Generates random matrices, distributed according to the G-Wishart distribution with parameters \eqn{b} and \eqn{D}, \eqn{W_G(b, D)} with respect to the graph structure \eqn{G}. 
Note this fuction works for both non-decomposable and decomposable graphs. 
}

\usage{ rgwish( n = 1, adj = NULL, b = 3, D = NULL, threshold = 1e-8 ) }

\arguments{
  \item{n}{number of samples required. }
  \item{adj}{
    adjacency matrix corresponding to the graph structure which can be non-decomposable or decomposable. It should be an upper triangular matrix in which \eqn{a_{ij}=1} 
    if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
    \code{adj} could be an object of class \code{"graph"}, from function \code{\link{graph.sim}}.
    It also could be an object of class "\code{sim}", from function \code{\link{bdgraph.sim}}.
    It also could be an object of class "\code{bdgraph}", from functions \code{\link{bdgraph.mpl}} or \code{\link{bdgraph}}.
    }
  \item{b}{degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. }
  \item{D}{positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. 
           The default is an identity matrix.}
  \item{threshold}{ threshold value for the convergence of sampling algorithm from G-Wishart.}
}

\details{
Sampling from G-Wishart distribution, \eqn{K \sim W_G(b, D)}, with density:

  \deqn{Pr(K) \propto |K| ^ {(b - 2) / 2} \exp \left\{- \frac{1}{2} \mbox{trace}(K \times D)\right\},}

which \eqn{b > 2} is the degree of freedom and \eqn{D} is a symmetric positive definite matrix.
}

\value{
A numeric array, say \eqn{A}, of dimension \eqn{(p \times p \times n)}, where each \eqn{A[,,i]}{A[,,i]} is a positive 
definite matrix, a realization of the G-Wishart distribution, \eqn{W_G(b, D)}.
Note, for the case \eqn{n=1}, the output is a matrix.
}

\references{
Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2:119-128, \doi{10.1002/sta4.23}

Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{gnorm}}, \code{\link{rwish}} }

\examples{
# Generating a 'circle' graph as a non-decomposable graph
adj <- graph.sim( p = 5, graph = "circle" )
adj    # adjacency of graph with 5 nodes
  
sample <- rgwish( n = 1, adj = adj, b = 3, D = diag( 5 ) )
round( sample, 2 ) 

sample <- rgwish( n = 5, adj = adj )
round( sample, 2 )  
}

\keyword{Wishart}
\keyword{ distribution }

