\name{summary.bdgraph}
\alias{summary.bdgraph}

\title{
 Summary function for "bdgraph" output
}
\description{
According to output of the BDMCMC algorithm, this function gives us the summary.  
}
\usage{
\method{summary}{bdgraph}(object, plot = TRUE, mode = "circle", edge.col = "black", label.col = "black", 
        vertex.col = "red", ...)
}

\arguments{
  \item{object}{
an object with S3 class "bdgraph". A list which is the result of the BDMCMC algorithm from the 'bdgraph' function. 
}
  \item{plot}{
logical: if TRUE (default) you will see the plot result.}
  \item{mode}{
the vertex placement algorithm which according to network package has three options "circle", "fruchtermanreingold", and "kamadakawai". 
The default is "circle".
}
  \item{edge.col}{
color for edges; may be given as a vector, adjacency matrix, or edge attribute
name, if edges are to be of different colors.
}
  \item{label.col}{
color for vertex labels; may be given as a vector or a vertex attribute name, if
labels are to be of different colors.
}
  \item{vertex.col}{
color for vertices; may be given as a vector or a vertex attribute name, if vertices
are to be of different colors.
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
  
  output <- bdgraph(data $ data, meanzero = T, iter = 2000)
  
  summary(output)
  
  summary(output, plot = FALSE)
  }
}
