\name{I_G}
\alias{I_G}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computing normalizing constant of G-Wishart distribution
}
\description{
Monte Carlo method for approximating the normalizing constant of G-Wishart
distribution.
}
\usage{
I_G(A, b, D, MC.iter = 300)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Upper truculer matrix which shows the starting graphs in which a_ij=1 if there 
is a link between notes i and j, otherwise a_ij=0  
}
  \item{b}{
Value for prior distribution of precision matrix
}
  \item{D}{
Positive definite matrix for prior distribution of precision matrix
}
  \item{MC.iter}{
Number of iterations for Monte Carlo approximation
}
}
\value{
Normalizing constant of G-Wishart distribution
}
\references{
Mohammadi, A. and E. Wit (2012). Efficient birth-death MCMC inference for 
Gaussian graphical models, Journal of the Royal Statistical Society: Series B,
submitted.

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable gaussian graphical models. Biometrika 
92(2), 317-335.
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\examples{
p=3
A=matrix(c(0,1,1,0,0,1,0,0,0),p,p,byrow=TRUE) # full graph
I_G(A, b=3, D=diag(p), MC.iter=50)
}

