\name{plot.simulate}
\alias{plot.simulate}

\title{ Plot function for \code{S3} class \code{"simulate"} }

\description{ Visualize structure of the true graph. }
\usage{
\method{plot}{simulate}( x, main = NULL, layout = layout.circle, ... )
}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"simulate"}, from function \code{\link{bdgraph.sim}}.  }
  
  \item{main}{ Graphical parameter (see plot).}
  
  \item{layout}{The vertex placement algorithm which is according to \pkg{igraph} package. The default is \code{"layout.circle"}.}
  
  \item{\dots}{System reserved (no specific usage).}
}
\references{
Mohammadi, A. and E. Wit (2014). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis} (Forthcoming). \url{http://arxiv.org/abs/1210.5371v6}

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108}

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling 1}, 225-230
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 10, p = 15 )

plot( data.sim )
}
}
