\name{I.g}
\alias{I.g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing normalizing constant of G-Wishart distribution
}
\description{
Monte Carlo method for approximating the normalizing constant of G-Wishart
distribution. The function uses the Monte Carlo method of Atay-Kayis and Massam (2005).
}
\usage{
I.g(A, b, D, MCiter = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
upper triangular matrix in which \eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}.}
  \item{b}{
the degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}{WG(b,D)}.
}
  \item{D}{
the positive definite matrix for G-Wishart distribution, \eqn{W_G(b,D)}{WG(b,D)}.
}
  \item{MCiter}{
the number of iterations for the Monte Carlo approximation.
}
}
\value{
the normalizing constant of G-Wishart distribution.
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, arXiv:1210.5371. \url{http://arxiv.org/abs/1210.5371}

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable gaussian graphical models. Biometrika 
92(2), 317-335.
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\examples{
A <- matrix(c(0,1,1,
              0,0,1,
		              0,0,0), 3, 3, byrow = TRUE) 
# matrix A shows full graph with 3 nodes and 3 edges
I.g(A, b=3, D=diag(3), MCiter=50)
}

