\name{print.simulate}
\alias{print.simulate}

\title{
 Print function for "bdgraph.sim" output
}
\description{
print the sample size, the number of nodes, number of linkes, type and sparcity of the true graph.  
}
\usage{
\method{print}{simulate}(x, ...)
}

\arguments{
  \item{x}{
an object with S3 class "simulate". A list which is the result from the "bdgraph.sim" function. 
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
\code{\link{bdgraph.sim}}
}
\examples{
\dontrun{
  data.sim <- bdgraph.sim(n = 10, p = 8, vis = T)
  
  print(data.sim)
  }
}
