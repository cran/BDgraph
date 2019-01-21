\name{gnorm}
\alias{gnorm}

\title{ Normalizing constant for G-Wishart }

\description{
Calculates log of the normalizing constant of G-Wishart distribution based on the Monte Carlo method, developed by Atay-Kayis and Massam (2005).
}

\usage{
gnorm( adj.g, b = 3, D = diag( ncol( adj.g ) ), iter = 100 )
}

\arguments{
  \item{adj.g}{
The adjacency matrix corresponding to the graph structure. It is an upper triangular matrix in which 
\eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, 
otherwise \eqn{a_{ij}=0}{aij = 0}.
}
  \item{b}{ The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. }
  
  \item{D}{
	The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, 
	\eqn{W_G(b,D)}. The default is an identity matrix.
	}
  \item{iter}{The number of iteration for the Monte Carlo approximation.}
}

\details{
Log of the normalizing constant approximation using Monte Carlo method for a G-Wishart distribution, \eqn{K \sim W_G(b, D)}, with density:

  \deqn{Pr(K) = \frac{1}{I(b, D)} |K| ^ {(b - 2) / 2} \exp \left\{- \frac{1}{2} \mbox{trace}(K \times D)\right\}.}
}

\value{ Log of the normalizing constant of G-Wishart distribution. }

\references{
Atay-Kayis, A. and Massam, H. (2005). A monte carlo method for computing the marginal likelihood in nondecomposable Gaussian graphical models, \emph{Biometrika}, 92(2):317-335

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Uhler, C., et al (2018) Exact formulas for the normalizing constants of Wishart distributions for graphical models, \emph{The Annals of Statistics} 46(1):90-118

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{ \code{\link{rgwish}}, \code{\link{rwish}} }

\examples{
\dontrun{
# adj.g: adjacency matrix of graph with 3 nodes and 2 links
adj.g <- matrix( c( 0, 0, 1,
                    0, 0, 1,
                    0, 0, 0 ), 3, 3, byrow = TRUE )		                
   
gnorm( adj.g, b = 3, D = diag( 3 ) )
}
}

