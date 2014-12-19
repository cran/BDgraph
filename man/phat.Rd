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
phat( output, round = 3 )
}

\arguments{
  \item{output}{An object with S3 class "bdgraph". }

  \item{round}{A number for rounding all probabilities to the specified number of decimal places (default 3).}
}
\value{
\item{phat}{Upper triangular matrix which shows the posterior inclusion probabilities for all possible edges.}
}
\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}} and \code{\link{bdgraph.sim}}}

\examples{
\dontrun{
# generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, graph = "circle", vis = TRUE )

output   <- bdgraph( data = data.sim )

phat( output, round = 2 )
}
}

