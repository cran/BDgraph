\name{traceplot}
\alias{traceplot}

\title{Trace plot of graph size}

\description{
Trace plot of graph size from all bdgraph output. It is a tool for checking the convergency of the BD-MCMC algorithm.
}

\usage{traceplot ( output, acf = FALSE, pacf = FALSE, main = NULL, ... )}

\arguments{
  \item{output}{ An object with S3 class "bdgraph". }
  \item{acf}{Logical: if TRUE you will see the plot of autocorrelation functions for graph sizes. The default is FALSE. }
  \item{pacf}{Logical: if TRUE you will see the plot of partial autocorrelations for graph sizes. The default is FALSE. }
  \item{main}{ Graphical parameter (see plot). }
  \item{\dots}{ System reserved (no specific usage). }
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian model selection in sparse Gaussian 
graphical models, arXiv:1210.5371v5. \url{http://arxiv.org/abs/1210.5371v5}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
  
  output   <- bdgraph( data = data.sim, iter = 2000, save.all = TRUE )
  
  traceplot(output)
  
  traceplot( output, acf = TRUE, pacf = TRUE )
  }
}
