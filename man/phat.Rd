\name{phat}
\alias{phat}

\title{ Posterior link probabilities }

\description{
Provides the posterior link probabilities for all possible links in the graph.
}

\usage{ phat( output, round = 3 ) }

\arguments{
  \item{output}{An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }

  \item{round}{A value for rounding all probabilities to the specified number of decimal places (default is 3).}
}

\value{
\item{phat}{Upper triangular matrix which corresponds the posterior probabilities for all possible links.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}} and \code{\link{bdgraph.sim}}}

\examples{
	\dontrun{
	# generating multivariate normal data from a 'circle' graph
	data.sim <- bdgraph.sim( n = 70, p = 6, graph = "circle", vis = TRUE )

	output   <- bdgraph( data = data.sim, iter = 10000 )

	phat( output, round = 2 )
	}
}

