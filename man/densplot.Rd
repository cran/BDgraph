\name{densplot}
\alias{densplot}

\title{
Plot of posterior distribution for graph size
}
\description{
This function is for ploting the posterior distribution of all possible graphs according to number
of their links.
}
\usage{
densplot(output, xlim = NULL, ylim = NULL, main = NULL)
}

\arguments{
  \item{output}{
a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{xlim}{
graphical parameter (see plot).
}
  \item{ylim}{
graphical parameter (see plot).
}
  \item{main}{
graphical parameter (see plot).
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
  
  densplot(output)
  }
}
