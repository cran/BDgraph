\name{rGWishart}
\alias{rGWishart}

\title{Sampling from G-Wishart distribution}

\description{
Generate n random matrices, distributed according to the G-Wishart distribution with parameters b and D, \eqn{W_G(b, D)}.
}

\usage{rGWishart(n = 1, b = 3, D = diag(2), G = NULL, method = "block gibbs", start.K = NULL)}

\arguments{
  \item{n}{integer sample size.}
  \item{b}{the degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}.}
  \item{D}{the positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}.}
  \item{G}{
adjacency matrix which shows the graph structure. It is an upper triangular matrix in which \eqn{g_{ij}=1} 
if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{g_{ij}=0}.
}
  \item{method}{
a character with two options "block gibbs" and "accept-reject" to determine a type of sampling method. Option 
"block gibbs" (default) determines a pair wise block gibbs sampling from G-Wishart distribution. Option 
"accept-reject" determines a accept-reject method for sampling from G-Wishart distribution.
}
  \item{start.K}{starting point for "block gibbs" method.}
}
\details{
Sampling from G-Wishart(b,D) distribution with density:

  \deqn{p(K) \propto |K| ^ {(b - 2) / 2} exp(- \frac{1}{2} trace(K \times D))}

which \eqn{b > 2} is the degree of freedom and D is a symmetric positive definite matrix.
}
\value{
a numeric array, say R, of dimension \eqn{(p \times p \times n)}, where each \eqn{R[,,i]}{R[,,i]} is a positive 
definite matrix, a realization of the G-Wishart distribution \eqn{W_G(b, D)}.
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}

Wang, H. and S.Z. Li (2012). Efficient Gaussian graphical model determination under G-Wishart prior distributions.
Electronic Journal of Statistics 6, 168-198.

Wang, H. and C. M. Carvalho (2010). Simulation of hyper-inverse Wishart distributions for non-decomposable graphs.
Electronic Journal of Statistics 4, 1470-1475.
}
\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
\dontrun{
   G <- toeplitz(c(0, 1, rep(0, 3)))
   G # graph with 5 nodes and 4 links
   
   sample <- rGWishart(n = 3, b = 3, D = diag(5), G = G, method = "block gibbs")
   round(sample, 2)
   }
}

