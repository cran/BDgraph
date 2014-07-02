\name{rgwish}
\alias{rgwish}

\title{Sampling from G-Wishart distribution}

\description{
Generate n random matrices, distributed according to the G-Wishart distribution with parameters b and D, \eqn{W_G(b, D)}.
}

\usage{rgwish( n = 1, G = NULL, b = 3, D = NULL )}

\arguments{
  \item{n}{The number of samples required. The default is 1.}
  \item{G}{
    Adjacency matrix corresponding to the graph structure. It is an upper triangular matrix in which \eqn{G_{ij}=1} 
    if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{G_{ij}=0}.
    }
  \item{b}{the degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. The default is 3.}
  \item{D}{the positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. 
           The default is an identity matrix.}
}
\details{
Sampling from G-Wishart(b,D) distribution with density:

  \deqn{p(K) \propto |K| ^ {(b - 2) / 2} exp(- \frac{1}{2} trace(K \times D))}

which \eqn{b > 2} is the degree of freedom and D is a symmetric positive definite matrix.
}
\value{
A numeric array, say A, of dimension \eqn{(p \times p \times n)}, where each \eqn{A[,,i]}{A[,,i]} is a positive 
definite matrix, a realization of the G-Wishart distribution \eqn{W_G(b, D)}.
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}

Lenkoski, A. (2013). A direct sampler for G-Wishart variates, Stat 2, 119-128.
}
\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
	\dontrun{
		G <- toeplitz( c( 0, 1, rep( 0, 3 ) ) )
		G # graph with 5 nodes and 4 links

		sample <- rgwish( n = 3, G = G, b = 3, D = diag(5) )
		round( sample, 2 )
	}
}

