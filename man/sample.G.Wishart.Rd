\name{sample.G.Wishart}
\alias{sample.G.Wishart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling from G-Wishart distribution
}
\description{
Tools for sampling from G-Wishart distribution according to Choleski 
decomposition of a Wishart variate with the identity as its scale parameter.
}
\usage{
sample.G.Wishart(A, b, D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
upper truculer matrix which show the starting graphs in which a_ij=1 if there 
is a link between notes i and j and otherwise a_ij=0  
}
  \item{b}{
value for prior distribution of precision matrix
}
  \item{D}{
positive definite matrix for prior distribution of precision matrix
}
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
p=5
b=3; D=diag(p)
A=matrix(0,p,p)
for (i in 1:(p-1)){
   for (j in (i+1):p){ 
      if (j==i+1) A[i,j]=1
      } 
   }
A[1,p]=1
round(sample.G.Wishart(A,b,D),2)
}
