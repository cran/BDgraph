\name{summary.bdgraph}
\alias{summary.bdgraph}

\title{ Summary function for \code{S3} class \code{"bdgraph"} }

\description{
This function provides a summary of the result from BD-MCMC sampling algorithm.  
}

\usage{
\method{summary}{bdgraph}( object, vis = TRUE, ... )
}

\arguments{
  \item{object}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{vis}{Logical: if TRUE (default) you will see the plot result.}
  \item{\dots}{System reserved (no specific usage).}
}

\value{
	\item{best.graph}{The adjacency matrix corresponding to the selected graph which has the highest posterior probability.}
	\item{phat}{Upper triangular matrix corresponding to the posterior probabilities for all possible links.}
	\item{Khat}{The estimated precision matrix.}
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling}, 1:225-230
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, vis = TRUE )

output   <- bdgraph( data = data.sim )

summary(output)

summary( output, vis = FALSE )
}
}
