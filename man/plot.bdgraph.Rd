\name{plot.bdgraph}
\alias{plot.bdgraph}

\title{
 Plot function for "bdgraph" output
}
\description{
According to output of the BD-MCMC algorithm, this function gives us the plot.  
}
\usage{
\method{plot}{bdgraph}(x, g = 1, layout = layout.circle, ...)
}

\arguments{
  \item{x}{
an object with S3 class "bdgraph". A list which is the result of the BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{g}{
you can select the number of graphical models with highest probabilities (default is 1).
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
  # generate the data (100 observations) from multivariate normal distribution
  data <- bdgraph.sim(n = 100, p = 8, vis = T)
  
  output <- bdgraph(data, meanzero = T, iter = 2000)
  
  plot(output)
  
  plot(output, g = 4)
  }
}
