\name{plot.sim}
\alias{plot.sim}

\title{
 Plot function for "bdgraph.sim" output
}
\description{
Visualize the structure of true graph.  
}
\usage{
\method{plot}{sim}(x, main = NULL, layout = layout.circle, ...)
}

\arguments{
  \item{x}{
an object with S3 class "sim". A list which is the result from the 'bdgraph.sim' function. 
}
  \item{main}{
graphical parameter(see plot).
}
  \item{layout}{
the vertex placement algorithm which is according to igraph package. 
The default is "layout.circle".
}
  \item{\dots}{
system reserved (no specific usege).
}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdgraph}}
}
\examples{
\dontrun{
  data <- bdgraph.sim(n = 10, p = 15, vis = T)
  
  plot(data)
  }
}
