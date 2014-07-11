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
I.g( G, b = 3, D = diag( ncol(G) ), mc = 100 )
}

\arguments{
  \item{G}{
Adjacency matrix corresponding to the graph structure. It is an upper triangular matrix in which 
\eqn{g_{ij}=1}{gij = 1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, 
otherwise \eqn{g_{ij}=0}{gij = 0}.
}
  \item{b}{ The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. The default is 3.}
  
  \item{D}{
The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, 
\eqn{W_G(b,D)}. The default is identity matrix.
}
  \item{mc}{The number of iteration for the Monte Carlo approximation. The default is 100.}
}
\details{
Normalizing constant approximation using Monte Carlo method for a G-Wishart(b,D):

  \deqn{p(K) = \frac{1}{I(b, D)} |K| ^ {(b - 2) / 2} exp(- \frac{1}{2} trace(K \times D))}
}
\value{
the normalizing constant of G-Wishart distribution.
}
\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}

Atay-Kayis, A. and H. Massam (2005). A monte carlo method for computing the 
marginal likelihood in nondecomposable Gaussian graphical models. Biometrika 
92(2), 317-335.
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\examples{
G <- matrix( c(0,0,1,
               0,0,1,
		               0,0,0), 3, 3, byrow = TRUE )
		                
# matrix G shows a graph with 3 nodes and 2 links
I.g( G, b = 3, D = diag(3) )
}

