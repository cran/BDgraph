\name{rcwish}
\alias{rcwish}

\title{ Sampling from Complex Wishart distribution }

\description{
Generates random matrices, distributed according to the complex Wishart distribution with parameters b and D, \eqn{CW(b, D)}.
}

\usage{ rcwish( n = 1, p = 2, b = 3, D = diag(p) ) }

\arguments{
  \item{n}{ The number of samples required. }
  \item{p}{ The number of variables (nodes). }
  \item{b}{ The degree of freedom for complex Wishart distribution, \eqn{CW(b, D)}. }
  \item{D}{ The positive definite \eqn{(p \times p)} "scale" matrix for complex Wishart distribution, \eqn{CW(b, D)}. 
            The default is an identity matrix. }
}

\details{
Sampling from complex Wishart distribution, \eqn{K \sim CW(b, D)}, with density:

  \deqn{Pr(K) \propto |K| ^ {b} \exp \left\{- \mbox{trace}(K \times D)\right\},}

which \eqn{b > 2} is the degree of freedom and D is a symmetric positive definite matrix.
}

\value{
A numeric array, say A, of dimension \eqn{(p \times p \times n)}, where each \eqn{A[,,i]}{A[,,i]} is a positive 
definite matrix, a realization of the complex Wishart distribution \eqn{CW(b, D)}.
}

\references{
Tank, A., N. Foti, and E. Fox (2015). Bayesian Structure Learning for Stationary Time Series, \emph{arXiv preprint arXiv:1505.03131}

D.R. Brillinger (2001). Time Series: Data Analysis and Theory, \emph{Holden-Day}

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 
}

\author{ Lang Liu, Nicholas Foti, Alex Tank and  Abdolreza Mohammadi }

\examples{
\dontrun{
p      <- 5
sample <- rcwish( n = 3, p = p, b = 3, D = diag(p) )
round( sample, 2 )  
}
}

