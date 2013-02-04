\name{I.g}
\alias{I.g}

\title{
Normalizing constant of G-Wishart distribution
}
\description{
Monte Carlo method for approximating the normalizing constant of G-Wishart
distribution. The function uses the Monte Carlo method of Atay-Kayis and Massam (2005).
}
\usage{
I.g(G, b = 3, D = diag(ncol(G)), MCiter = 100)
}

\arguments{
  \item{G}{
adjacency matrix which shows the garph structure. It is an upper triangular matrix in which 
\eqn{a_{ij}=1}{Gij = 1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, 
otherwise \eqn{G_{ij}=0}{Gij = 0}.
}
  \item{b}{
the degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}{WG(b, D)}. The default is 3.
}
  \item{D}{
the positive definite \eqn{(p * p)}{(p * p)} "scale" matrix for G-Wishart distribution, 
\eqn{W_G(b,D)}{WG(b,D)}. The default is identity matrix.
}
  \item{MCiter}{
the number of iterations for the Monte Carlo approximation. The default is 100.
}
}
\details{
normalizing constant approximation using Monte Carlo method for a G-Wishart(b,D):

  \deqn{p(K) = \frac{1}{I(b, D)} |K| ^ {(b - 2) / 2} exp(- \frac{1}{2} trace(K \times D))}
}
\value{
the normalizing constant of G-Wishart distribution.
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable gaussian graphical models. Biometrika 
92(2), 317-335.
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\examples{
G <- matrix(c(0,0,1,
              0,0,1,
		              0,0,0), 3, 3, byrow = TRUE) 
# matrix G shows full graph with 3 nodes and 3 edges
I.g(G, b = 3, D = diag(3))
}

