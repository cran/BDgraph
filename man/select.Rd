\name{select}
\alias{select}

\title{ Graph selection }

\description{
Provides the selected graph which, based on input, could be a graph with links for which their estimated posterior probabilities are greater than 0.5 (as a default)
or a graph with the highest posterior probability; see examples.
%For more specific selection of graphs consult the 'prob' function.  
}

\usage{ select( bdgraph.obj, cut = NULL, vis = FALSE ) }

\arguments{
  \item{bdgraph.obj}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  \item{cut}{ Threshold for including the links in the selected graph based on the estimated posterior probabilities of the links; see the examples. }  
  \item{vis}{ Visualize the selected graph structure. The default value is FALSE. }
}

\value{	\item{G}{ An adjacency matrix corresponding to the selected graph.} }

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Gaussian Copula Graphical Modeling for Dupuytren Disease, \emph{arXiv:1501.04849} 
}

\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
  
bdgraph.obj <- bdgraph( data = data.sim )
   
select( bdgraph.obj )
  
bdgraph.obj <- bdgraph( data = data.sim, save.all = TRUE )
  
select( bdgraph.obj )
  
select( bdgraph.obj, cut = 0.5, vis = TRUE )
}
}
