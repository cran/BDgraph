\name{prob}
\alias{prob}

\title{Posterior probabilities of the graphs}

\description{
According to the output of the BD-MCMC algorithm, this function gives the 
posterior probabilities for the 'g' most likely graphs or a specific graph, 'G'.
}

\usage{prob( output, g = 4, G = NULL )}

\arguments{
  \item{output}{ An object with S3 class "bdgraph". }
  \item{g}{ The number of graphs with the highest probabilities to be shown (default is 4).
            This option is ignored if 'G' is specified. }
  \item{G}{ Adjacency matrix corresponding to a graph structure. It is an upper triangular matrix in which 
            \eqn{G_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{G_{ij}=0}.
            It also can be an object with S3 class "simulate". }
}
\value{
\item{ best.G }{ The graphs with the highest posterior probabilities.}

\item{ prob.G }{ A vector which includes posterior probabilities of the graphs in 'best.G'.}
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{ \code{\link{bdgraph}} }

\examples{
\dontrun{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, size = 6, vis = TRUE )
  
  output   <- bdgraph( data = data.sim, iter = 2000 )
 
  # Estimated posterior probability of the true graph
  prob( output, G = data.sim )
  
  # Estimated posterior probability of the first and second graphs with highest probabilities
  prob( output, g = 2 )
  }
}
