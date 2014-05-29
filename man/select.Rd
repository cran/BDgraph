\name{select}
\alias{select}

\title{Selecting the best graph}

\description{
According to output of the BD-MCMC algorithm, this function gives the graphs with the highest posterior probabilities.
For more specific selection of graphs consult the 'prob' function.  
}

\usage{select( output, vis = FALSE )}

\arguments{
  \item{output}{ An object with S3 class "bdgraph". }
  
  \item{vis}{Logical: if TRUE you will see the plot of best graph. The default is FALSE. }
}
\value{
\item{G}{ Adjacency matrix corresponding to the best graph (graph with the highest posterior probability).}
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
  
  output   <- bdgraph( data = data.sim, iter = 2000 )
  
  select(output)
  }
}
