\name{plot.bdgraph}
\alias{plot.bdgraph}

\title{
 Plot function for "bdgraph" output
}
\description{
According to the output of our BD-MCMC algorithm, this function gives us the plot.  
}
\usage{
\method{plot}{bdgraph}( x, g = 1, layout = layout.circle, ... )
}

\arguments{
  \item{x}{ An object with S3 class "bdgraph". }
  \item{g}{The number of graphical models with highest probabilities (default is 1).}
  
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
\dontrun{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, vis = TRUE )
  
  output   <- bdgraph( data = data.sim, iter = 2000 )
  
  plot(output)
  
  plot( output, g = 4 )
  }
}
