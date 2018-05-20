\name{plotroc}
\alias{plotroc}

\title{ ROC plot }

\description{
Draws the ROC curve according to the true graph structure for object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}.
}

\usage{ 
plotroc( sim.obj, bdgraph.obj, bdgraph.obj2 = NULL, bdgraph.obj3 = NULL, 
         bdgraph.obj4 = NULL, cut = 20, smooth = FALSE, label = TRUE, 
         main = "ROC Curve" ) 
}

\arguments{
  \item{sim.obj}{
	An object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}.
	It also can be the adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and
	\eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
  }
  
  \item{bdgraph.obj}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{ssgraph} of \code{R} package \pkg{ssgraph}.
					  It also can be an object of \code{S3} class \code{"select"}, from the function \code{huge.select} of \code{R} package \pkg{huge}.
                      It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links. }
  
  \item{bdgraph.obj2}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{ssgraph} of \code{R} package \pkg{ssgraph}.
					  It also can be an object of \code{S3} class \code{"select"}, from the function \code{huge.select} of \code{R} package \pkg{huge}.
                       It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
                       It is for comparing two different approaches. }
  
  \item{bdgraph.obj3}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{ssgraph} of \code{R} package \pkg{ssgraph}.
					  It also can be an object of \code{S3} class \code{"select"}, from the function \code{huge.select} of \code{R} package \pkg{huge}.
                       It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
                       It is for comparing three different approaches. }
  
  \item{bdgraph.obj4}{ An object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. 
					  It also can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{ssgraph} of \code{R} package \pkg{ssgraph}.
					  It also can be an object of \code{S3} class \code{"select"}, from the function \code{huge.select} of \code{R} package \pkg{huge}.
                       It also can be an upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
                       It is for comparing four different approaches. }
  
  \item{cut}{ Number of cut points.}
  
  \item{smooth}{ Logical: for smoothing the ROC curve.} 
  \item{label}{ Logical: for adding legend to the ROC plot.} 
  \item{main}{ An overall title for the plot.} 
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{arXiv preprint arXiv:1501.05108} 

Dobra, A. and A. Mohammadi (2017). Loglinear Model Selection and Human Mobility, \emph{arXiv preprint arXiv:1711.02623}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C} 

Mohammadi, A., Massam H., and G. Letac (2017). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416} 
}

\author{ Reza Mohammadi and Ernst Wit }

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
legend( "bottomright", c( "GGMs", "GCGMs" ), lty = c( 1,2 ), col = c( "black", "red" ) )   
}
}
