\name{print.simulate}
\alias{print.simulate}

\title{ Print function for \code{S3} class \code{"simulate"} }

\description{ Print the information about the type of data, the sample size, the graph type, the number of nodes, number of linkes and sparsity of the true graph. }

\usage{ \method{print}{simulate}( x, ... ) }

\arguments{
  \item{x}{ An object of \code{S3} class \code{"simulate"}, from function \code{\link{bdgraph.sim}}.  }
  \item{\dots}{ System reserved (no specific usege). }
}
\references{
Mohammadi, A. and E. Wit (2015). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling}, 1:225-230
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdgraph.sim}}
}
\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph 
data.sim <- bdgraph.sim( n = 20, p = 10, vis = TRUE )

print( data.sim )
}
}
