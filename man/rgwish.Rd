\name{rgwish}
\alias{rgwish}

\title{Sampling from G-Wishart distribution}

\description{
Generates random matrices, distributed according to the G-Wishart distribution with parameters b and D, \eqn{W_G(b, D)} with respect to the graph structure \eqn{G}. 
Note this fuction works for both non-decomposable and decomposable graphs. 
}

\usage{ rgwish( n = 1, adj = NULL, b = 3, D = NULL, threshold = 1e-8 ) }

\arguments{
  \item{n}{The number of samples required. }
  \item{adj}{
    The adjacency matrix corresponding to the graph structure which can be non-decomposable or decomposable. It should be an upper triangular matrix in which \eqn{a_{ij}=1} 
    if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
    \code{adj} could be an object of class \code{"graph"}, from function \code{\link{graph.sim}}.
    It also could be an object of class \code{"sim"}, from function \code{\link{bdgraph.sim}}.
    It also could be an object of class \code{"bdgraph"}, from functions \code{\link{bdgraph.mpl}} or \code{\link{bdgraph}}.
    }
  \item{b}{The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. }
  \item{D}{The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. 
           The default is an identity matrix.}
  \item{threshold}{ The threshold value for the convergence of sampling algorithm from G-Wishart.}
}

\details{
Sampling from G-Wishart distribution, \eqn{K \sim W_G(b, D)}, with density:

  \deqn{Pr(K) \propto |K| ^ {(b - 2) / 2} \exp \left\{- \frac{1}{2} \mbox{trace}(K \times D)\right\},}

which \eqn{b > 2} is the degree of freedom and D is a symmetric positive definite matrix.
}

\value{
A numeric array, say A, of dimension \eqn{(p \times p \times n)}, where each \eqn{A[,,i]}{A[,,i]} is a positive 
definite matrix, a realization of the G-Wishart distribution, \eqn{W_G(b, D)}.
Note, for the case \eqn{n=1}, the output is a matrix.
}

\references{
Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2:119-128

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and Wit, E. C. (2017). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108v5} 

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{gnorm}}, \code{\link{rwish}}, \code{\link{rgcwish}} }

\examples{
# Generating a 'circle' graph as a non-decomposable graph
adj <- graph.sim( p = 5, graph = "circle" )
adj    # adjacency of graph with 5 nodes
  
sample <- rgwish( n = 1, adj = adj, b = 3, D = diag( 5 ) )
round( sample, 2 ) 

sample <- rgwish( n = 5, adj = adj )
round( sample, 2 )  
}

