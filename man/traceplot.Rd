\name{traceplot}
\alias{traceplot}

\title{
Trace plot of graph size from bdgraph output
}
\description{
Trace plot of graph size from all bdgraph output. 
It is a tool for checking the convergency of the BDMCMC algorithm.
}
\usage{
traceplot (output, acf = FALSE, pacf = FALSE, xlim = NULL, ylim = NULL, main = NULL)
}

\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdgraph' function. 
}
  \item{acf}{
logical: if TRUE you will see the plot of autocorrelation functions for graph sizes. The default is FALSE. 
}
  \item{pacf}{
logical: if TRUE you will see the plot of partial autocorrelations for graph sizes. The default is FALSE. 
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
  
  output <- bdgraph(data $ data, meanzero = T, iter = 2000, save.all = TRUE)
  
  traceplot(output)
  
  traceplot(output, acf = T, pacf = T)
  }
}
