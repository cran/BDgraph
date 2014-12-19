\name{rwish}
\alias{rwish}

\title{Sampling from Wishart distribution}

\description{
Generate n random matrices, distributed according to the Wishart distribution with parameters b and D, \eqn{W(b, D)}.
}

\usage{rwish( n = 1, p = 2, b = 3, D = diag(p) )}

\arguments{
  \item{n}{ The number of samples required. The default is 1. }
  \item{p}{ The number of variables (nodes). The default value is 2. }
  \item{b}{ The degree of freedom for Wishart distribution, \eqn{W(b, D)}. The default is 3. }
  \item{D}{ The positive definite \eqn{(p \times p)} "scale" matrix for Wishart distribution, \eqn{W(b, D)}. 
            The default is an identity matrix. }
}
\details{
Sampling from Wishart(b,D) distribution with density:

  \deqn{p(K) \propto |K| ^ {(b - 2) / 2} exp(- \frac{1}{2} trace(K \times D))}

which \eqn{b > 2} is the degree of freedom and D is a symmetric positive definite matrix.
}
\value{
A numeric array, say A, of dimension \eqn{(p \times p \times n)}, where each \eqn{A[,,i]}{A[,,i]} is a positive 
definite matrix, a realization of the Wishart distribution \eqn{W(b, D)}.
}
\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
\dontrun{
sample <- rwish( n = 3, p = 5, b = 3, D = diag(p) )
sample
}
}

