\name{prob}
\alias{prob}

\title{Posterior probabilities of the graphs}

\description{
According to the output of the BD-MCMC algorithm, this function gives us the 
posterior probability for the most highest graphs or a specific graph.
}

\usage{prob(output, g = 2, G = NULL)}

\arguments{
  \item{output}{a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. }
  \item{g}{you can select the number of graphs with highest probabilities (default is 1). }
  \item{G}{ adjacency matrix which shows the graph structure. It is an upper triangular matrix in which 
  \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}. }
}
\value{
\item{ best.G }{ the graphs with the highest posterior probabilities.}

\item{ prob.G }{ a vector which includes posterior probabilities of the graphs in 'best.G'.}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{ \code{\link{bdgraph}} }

\examples{
\dontrun{
  # generating the data from a random graph
  data.sim <- bdgraph.sim(n = 100, p = 8, vis = TRUE)
  
  output   <- bdgraph(data.sim, mean = 0, iter = 2000)
  
  prob(output, G = data.sim)
  
  prob(output, g = 2)
  }
}
