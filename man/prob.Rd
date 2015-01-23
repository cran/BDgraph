\name{prob}
\alias{prob}

\title{Posterior probabilities of the graphs}

\description{
Provides the posterior probabilities for the \code{'g'} most likely graphs or a specific graph, \code{'G'}.
}

\usage{prob( output, g = 4, G = NULL )}

\arguments{
  \item{output}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{g}{ The number of graphs with the highest posterior probabilities to be shown (default is 4).
            This option is ignored if \code{'G'} is specified. }
  \item{G}{ Adjacency matrix corresponding to a graph structure. It is an upper triangular matrix in which 
            \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}.
            It also can be an object of \code{S3} class \code{"simulate"}, from function \code{\link{bdgraph.sim}}. }
}
\value{
\item{ best.G }{ The graphs with the highest posterior probabilities.}

\item{ prob.G }{ A vector of the posterior probabilities of the graphs corresponding to \code{'best.G'}.}
}
\references{
Mohammadi, A. and E. Wit (2014). Bayesian structure learning in sparse Gaussian graphical models, \emph{Bayesian Analysis} (Forthcoming). \url{http://arxiv.org/abs/1210.5371v6}

Mohammadi, A. and E. Wit (2015). BDgraph: Bayesian Undirected Graph Estimation in R, \emph{Arxiv preprint arXiv:1501.05108}

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Copula Gaussian Graphical Models, \emph{Arxiv preprint arXiv:1501.04849} 

Mohammadi, A., F. Abegaz Yazew, and E. Wit (2014). Bayesian Copula Gaussian graphical modelling, \emph{Proceedings of the 29th International Workshop on Statistical Modelling 1}, 225-230
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{ \code{\link{bdgraph}} }

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 6, vis = TRUE )

output   <- bdgraph( data = data.sim )

# Estimated posterior probability of the true graph
prob( output, G = data.sim )

# Estimated posterior probability of the first and second graphs with highest probabilities
prob( output, g = 2 )
}
}
