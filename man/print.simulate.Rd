\name{print.simulate}
\alias{print.simulate}

\title{
 Print function for "bdgraph.sim" output
}
\description{
Print the sample size, the number of nodes, number of linkes, type and sparsity of the true graph.  
}
\usage{
\method{print}{simulate}( x, ... )
}

\arguments{
  \item{x}{ An object with S3 class "simulate". }
  \item{\dots}{ System reserved (no specific usege). }
}
\references{
Mohammadi, A. and E. C. Wit (2013). Bayesian structure learning in sparse Gaussian 
graphical models, arXiv:1210.5371v6. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdgraph.sim}}
}
\examples{
	\dontrun{
		# generating synthetic multivariate normal data from a 'random' graph 
		data.sim <- bdgraph.sim( n = 10, p = 8, vis = TRUE )

		print(data.sim)
	}
}
