\name{rwish}
\alias{rwish}

\title{ Sampling from Wishart distribution }

\description{
Generates random matrices, distributed according to the Wishart distribution with parameters \eqn{b} and \eqn{D}, \eqn{W(b, D)}.
}

\usage{ rwish( n = 1, p = 2, b = 3, D = diag( p ) ) }

\arguments{
  \item{n}{ The number of samples required. }
  \item{p}{ The number of variables (nodes). }
  \item{b}{ The degree of freedom for Wishart distribution, \eqn{W(b, D)}. }
  \item{D}{ The positive definite \eqn{(p \times p)} "scale" matrix for Wishart distribution, \eqn{W(b, D)}. 
            The default is an identity matrix. }
}

\details{
Sampling from Wishart distribution, \eqn{K \sim W(b, D)}, with density:

  \deqn{Pr(K) \propto |K| ^ {(b - 2) / 2} \exp \left\{- \frac{1}{2} \mbox{trace}(K \times D)\right\},}

which \eqn{b > 2} is the degree of freedom and \eqn{D} is a symmetric positive definite matrix.
}

\value{
A numeric array, say A, of dimension \eqn{(p \times p \times n)}, where each \eqn{A[,,i]}{A[,,i]} is a positive 
definite matrix, a realization of the Wishart distribution \eqn{W(b, D)}.
Note, for the case \eqn{n=1}, the output is a matrix.
}

\references{
Lenkoski, A. (2013). A direct sampler for G-Wishart variates, \emph{Stat}, 2:119-128

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{gnorm}}, \code{\link{rgwish}}, \code{\link{rgcwish}} }

\examples{
sample <- rwish( n = 3, p = 5, b = 3, D = diag( 5 ) )

round( sample, 2 )  
}

\keyword{Wishart}
\keyword{distribution}
