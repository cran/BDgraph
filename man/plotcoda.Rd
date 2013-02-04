\name{plotcoda}
\alias{plotcoda}

\title{Convergency plots}

\description{
This function according to "bdgraph" output gives us several plots. These plots are  
useful tools for checking the convergency of the BD-MCMC algorithm.
}

\usage{plotcoda(output, thin = NULL, trace = TRUE, main = NULL, ...)}

\arguments{
  \item{output}{a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. }
  
  \item{thin}{it is a option for getting fast result for a cumulative plot according to part of the iterations.}
  
  \item{trace}{logical: if TRUE (default), tracing information is printed.}
  
  \item{main}{ graphical parameter (see plot). }
  
  \item{\dots}{ system reserved (no specific usage). }  
}
\details{
Note that a spending time for this function depends on the graph. It should be slow for the 
high-dimensional graphs. To make it faster you can choose bigger value for 'thin'.
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
  data.sim <- bdgraph.sim(n = 100, p = 8, graph = "circle", vis = TRUE)
  
  output <- bdgraph(data.sim, mean = 0, iter = 2000, save.all = TRUE)
 
  plotcoda(output)
  }
}
