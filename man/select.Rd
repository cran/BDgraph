\name{select}
\alias{select}

\title{Selecting the best graph}

\description{
Provides the graphs with the highest posterior probabilities.
For more specific selection of graphs consult the 'prob' function.  
}

\usage{ select( output, cut = NULL, vis = FALSE ) }

\arguments{
  \item{output}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{cut}{ Thrashold for including the links in the selected graph based on the posterior probabilities of the links; See example. }  
  \item{vis}{Logical: if TRUE you will see the plot of best graph. The default is FALSE. }
}

\value{
	\item{G}{ Adjacency matrix corresponding to the selected graph.}
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
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

output   <- bdgraph( data = data.sim )

select( output )

select( output, cut = 0.5, vis = TRUE )
}
}
