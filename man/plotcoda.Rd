\name{plotcoda}
\alias{plotcoda}

\title{
Cumulative occupancy fractions for checking the convergency
}
\description{
Plot the cumulative occupancy fractions for all possible edges. 
It is a tool for checking the convergency of the BDMCMC algorithm.
}
\usage{
plotcoda(output, skip = ceiling(length(output $ allA) / 2000), verbose = TRUE, main = NULL)
}

\arguments{
  \item{output}{
a list which is the result of BDMCMC algorithm from the 'bdgraph' function. 
}
  \item{skip}{
it is a option for getting fast result for a cumulative plot according to part of the iterations.
}
  \item{verbose}{
logical: if verbose = TRUE, tracing information printing is abled. The default value is TRUE.
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
  # generate the data (100 observations) from multivariate normal distribution
  data <- bdgraph.sim(n = 100, p = 8, vis = T)
  
  output <- bdgraph(data $ data, meanzero = T, iter = 2000, save.all = TRUE)
 
  plotcoda(output)
  }
}
