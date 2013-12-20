\name{traceplot}
\alias{traceplot}

\title{Trace plot of graph size}

\description{
Trace plot of graph size from all bdgraph output. It is a tool for checking the convergency of the BD-MCMC algorithm.
}

\usage{traceplot (output, acf = FALSE, pacf = FALSE, main = NULL, ...)}

\arguments{
  \item{output}{a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. }
  \item{acf}{logical: if TRUE you will see the plot of autocorrelation functions for graph sizes. The default is FALSE. }
  \item{pacf}{logical: if TRUE you will see the plot of partial autocorrelations for graph sizes. The default is FALSE. }
  \item{main}{ graphical parameter (see plot). }
  \item{\dots}{ system reserved (no specific usage). }
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
  
  output   <- bdgraph(data.sim, mean = 0, iter = 2000, save.all = TRUE)
  
  traceplot(output)
  
  traceplot(output, acf = TRUE, pacf = TRUE)
  }
}
