\name{plot.bdgraph}
\alias{plot.bdgraph}

\title{ Plot function for \code{S3} class \code{"bdgraph"} }

\description{ Visualize pattern of the graphs with the highest posterior probabilities. }

\usage{
\method{plot}{bdgraph}( x, g = 1, layout = layout.circle, ... )
}

\arguments{
  \item{x}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{g}{ The number of graphs with highest probabilities (default is 1).}
  
  \item{layout}{The vertex placement algorithm which is according to igraph package. The default is \code{"layout.circle"}.}
  
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
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, vis = TRUE )

output   <- bdgraph( data = data.sim )

plot(output)

plot( output, g = 4 )
}
}
