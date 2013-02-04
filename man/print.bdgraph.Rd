\name{print.bdgraph}
\alias{print.bdgraph}

\title{
 Print function for "bdgraph" output
}
\description{
According to output of the BD-MCMC algorithm, this function gives us the print of the result.  
}
\usage{
\method{print}{bdgraph}(x, round = 3, Khat = FALSE, phat = FALSE, ...)
}

\arguments{
  \item{x}{
an object with S3 class "bdgraph". A list which is the result of the BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{round}{
a number for rounding all probabilities to the specified number of decimal places (default 3).
}
  \item{Khat}{
logical: if TRUE you will see the estimation of precisiom matrix. The default is FALSE.
}
  \item{phat}{
logical: if TRUE you will see the posterior edge inclusion probabilities. The default is FALSE.
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
  
  print(output)
  
  print(output, Khat = TRUE, phat = TRUE)
  }
}
