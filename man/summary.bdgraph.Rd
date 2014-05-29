\name{summary.bdgraph}
\alias{summary.bdgraph}

\title{Summary function for "bdgraph" output}

\description{
This function gives us a summary of the result from BD-MCMC algorithm according to output from "bdgraph" function.  
}
\usage{
\method{summary}{bdgraph}( object, vis = TRUE, layout = layout.circle, ... )
}

\arguments{
  \item{object}{ An object with S3 class "bdgraph". }
  \item{vis}{Logical: if TRUE (default) you will see the plot result.}
  \item{layout}{ The vertex placement algorithm which is according to igraph package. The default is "layout.circle". }
  \item{\dots}{ System reserved (no specific usage).}
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
  
  summary(output)
  
  summary( output, vis = FALSE )
  }
}
