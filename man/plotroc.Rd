\name{plotroc}
\alias{plotroc}

\title{ ROC plot }

\description{
Draws the ROC curve according to the true graph structure for object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}.
}

\usage{ plotroc( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL,
                 cut.num = 20, smooth = FALSE, label = TRUE ) }

\arguments{
  \item{sim.obj}{
	An object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}.
	It also can be the adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and
	\eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
  }
  
  \item{bdgraph.obj}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
                      It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links. }
  
  \item{bdgraph.obj2}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
                       It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
                       It is for comparing two different approaches. }
  
  \item{bdgraph.obj3}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
                       It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
                       It is for comparing three different approaches. }
  
  \item{cut.num}{ Number of cut points. The default value is 20. }
  
  \item{smooth}{ Logical: for smoothing the ROC curve. The default is FALSE.} 
  \item{label}{ Logical: for adding legend to the ROC plot. The default is TRUE.} 
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv:1501.05108} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2016). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{\code{\link{bdgraph}} and \code{\link{compare}}}

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 30, p = 6, size = 7, vis = TRUE )
   
# Runing sampling algorithm
bdgraph.obj <- bdgraph( data = data.sim, iter = 10000 )
# Comparing the results
plotroc( data.sim, bdgraph.obj )
   
# To compare the results based on CGGMs approach
bdgraph.obj2 <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )
# Comparing the resultss
plotroc( data.sim, bdgraph.obj, bdgraph.obj2, label = FALSE )
legend( "bottomright", c( "GGMs", "GCGMs" ), lty = c( 1,2 ), col = c( 1, 4 ) )   
}
}
