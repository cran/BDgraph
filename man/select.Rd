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
  \item{cut}{ Threshold for including the links in the selected graph based on the posterior probabilities of the links; See the examples. }  
  \item{vis}{Logical: if TRUE you will see the plot of best graph. The default is FALSE. }
}

\value{	\item{G}{ Adjacency matrix corresponding to the selected graph.} }

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 
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
