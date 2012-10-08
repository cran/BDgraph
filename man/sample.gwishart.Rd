\name{sample.gwishart}
\alias{sample.gwishart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling from G-Wishart distribution
}
\description{
Tools for sampling from G-Wishart distribution according to Choleski 
decomposition of a Wishart variate with the identity as its scale parameter.
The function uses the sampling method of Atay-Kayis and Massam (2005).
}
\usage{
sample.gwishart(A, b, D, round = 3)
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
  \item{round}{
a number for rounding all values to the specified number of decimal places (default 3).
}
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
p <- 5
A <- matrix(0, p, p)
for (i in 1:(p-1)){
   for (j in (i+1):p){ 
      if (j==i+1) A[i,j] <- 1
      } 
   }
A[1,p] <- 1
# matrix A shows circle graph with 5 nodes and 5 edges
sample.gwishart(A, b=3, D=diag(p), round=1)
}
