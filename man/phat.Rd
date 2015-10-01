\name{phat}
\alias{phat}

\title{ Estimated posterior link probabilities }

\description{ 
Provides the estimated posterior link probabilities for all possible links in the graph.
}

\usage{ phat( x, round = 3 ) }

\arguments{
  \item{x}{An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }

  \item{round}{A value for rounding all probabilities to the specified number of decimal places (default is 3).}
}

\value{
\item{phat}{ An upper triangular matrix which corresponds the estimated posterior probabilities for all possible links.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}} and \code{\link{bdgraph.sim}}}

\examples{
	\dontrun{
	# Generating multivariate normal data from a 'circle' graph
	data.sim <- bdgraph.sim( n = 70, p = 6, graph = "circle", vis = TRUE )
   
	output   <- bdgraph( data = data.sim, iter = 10000 )
  
	phat( output, round = 2 )
	}
}

