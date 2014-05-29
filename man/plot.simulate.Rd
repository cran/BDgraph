\name{plot.simulate}
\alias{plot.simulate}

\title{
 Plot function for "bdgraph.sim" output
}
\description{
Visualize the structure of true graph.  
}
\usage{
\method{plot}{simulate}( x, main = NULL, layout = layout.circle, ... )
}

\arguments{
  \item{x}{An object with S3 class "simulate".  }
  
  \item{main}{ Graphical parameter (see plot).}
  
  \item{layout}{The vertex placement algorithm which is according to igraph package. The default is "layout.circle".}
  
  \item{\dots}{System reserved (no specific usage).}
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 10, p = 15 )
  
  plot(data.sim)
}
