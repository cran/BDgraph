\name{plotroc}
\alias{plotroc}

\title{ ROC plot }

\description{
Draws the receiver operating characteristic (ROC) curve according to the true graph structure for object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}.
}

\usage{ 
plotroc( target, est, est2 = NULL, est3 = NULL, est4 = NULL, 
         cut = 20, smooth = FALSE, label = TRUE, main = "ROC Curve" ) 
}

\arguments{
  \item{target}{
    adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
  }
  
  \item{est,
        est2,
        est3,
        est4}{
        upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
        It can be an object with \code{S3} class "\code{bdgraph}" from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be an object of \code{S3} class \code{"select"}, from the function \code{\link[huge]{huge.select}} of \code{R} package \code{\link[huge]{huge}}.
        Options \code{est2, est3} and \code{est4} are for comparing two or more different approaches.
        }

  \item{cut   }{ number of cut points.}
  \item{smooth}{ logical: for smoothing the ROC curve.} 
  \item{label }{ logical: for adding legend to the ROC plot.} 
  \item{main  }{ overall title for the plot.} 
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
  
Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} }

\seealso{\code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{compare}} }

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

legend( "bottomright", c( "GGMs", "GCGMs" ), lty = c( 1, 2 ), col = c( "black", "red" ) )   
}
}

\keyword{hplot}

