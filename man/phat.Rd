\name{phat}
\alias{phat}

\title{
Posterior edge inclusion probabilities
}
\description{
According to the output of BD-MCMC algorithm, this function gives us the posterior edge inclusion 
probabilities for all possible edges.
}
\usage{
phat(output, round = 3)
}

\arguments{
  \item{output}{
a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{round}{
a number for rounding all probabilities to the specified number of decimal places (default 3).
}
}
\value{
\item{phat}{upper triangular matrix which shows the posterior inclusion probabilities for all possible edges.}
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
  data <- bdgraph.sim(n = 100, p = 8, graph = "circle", vis = T)
  
  output <- bdgraph(data, meanzero = T, iter = 2000)
  
  phat(output, round = 2)
  }
}

