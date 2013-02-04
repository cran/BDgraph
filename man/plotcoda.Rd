\name{plotcoda}
\alias{plotcoda}

\title{
Convergency plots
}
\description{
This funciton according to "bdgraph" output gives us several plots. These plots are  
useful tools for checking the convergency of the BD-MCMC algorithm.
}
\usage{
plotcoda(output, skip = NULL, trace = TRUE, xlim = NULL, ylim = NULL, main = NULL)
}

\arguments{
  \item{output}{
a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{skip}{
it is a option for getting fast result for a cumulative plot according to part of the iterations.
}
  \item{trace}{
logical: if trace = TRUE, tracing information printing is abled. The default value is TRUE.
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
\details{
Note that a spending time for this function depences on the graph. It should be slow for the 
high-dimensional graphs. To make it faster you can choose bigger value for 'skip'.
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
  # generate the data from multivariate normal distribution
  data.sim <- bdgraph.sim(n = 100, p = 8, graph = "circle", vis = T)
  
  output <- bdgraph(data.sim, meanzero = T, iter = 2000, save.all = TRUE)
 
  plotcoda(output)
  }
}
