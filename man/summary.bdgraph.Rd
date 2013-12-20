\name{summary.bdgraph}
\alias{summary.bdgraph}

\title{Summary function for "bdgraph" output}

\description{
This function gives us a summary of the result from BD-MCMC algorithm according to output from "bdgraph" function.  
}
\usage{
\method{summary}{bdgraph}(object, vis = TRUE, layout = layout.circle, ...)
}

\arguments{
  \item{object}{
an object with S3 class "bdgraph". A list which is the result of the BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{vis}{logical: if TRUE (default) you will see the plot result.}
  \item{layout}{
the vertex placement algorithm which is according to igraph package. The default is "layout.circle".
}
  \item{\dots}{ system reserved (no specific usage).}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
  # generating the data from a random graph
  data.sim <- bdgraph.sim(n = 100, p = 8, vis = TRUE)
  
  output   <- bdgraph(data.sim, mean = 0, iter = 2000)
  
  summary(output)
  
  summary(output, vis = FALSE)
  }
}
