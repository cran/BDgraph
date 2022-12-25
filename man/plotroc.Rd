\name{plotroc}
\alias{plotroc}

\title{ ROC plot }

\description{
Draws the receiver operating characteristic (ROC) curve according to the true graph structure for object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}.
}

\usage{ 
plotroc( pred, actual, cut = 200, smooth = FALSE, calibrate = TRUE, 
         linetype = NULL, color = NULL, size = 1, main = "ROC Curve", 
         xlab = "False Postive Rate", ylab = "True Postive Rate",
         legend = TRUE, legend.size = 17, legend.position = c( 0.7, 0.3 ),
         labels = NULL, auc = TRUE, theme = ggplot2::theme_minimal() ) 
}

\arguments{

  \item{pred}{
        upper triangular matrix corresponding to the estimated posterior probabilities for all possible links.
        It can be an object with \code{S3} class "\code{bdgraph}" from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be an object of \code{S3} class \code{"select"}, from the function \code{\link[huge]{huge.select}} of \code{R} package \code{\link[huge]{huge}}.
        Options \code{est2, est3} and \code{est4} are for comparing two or more different approaches.
        }

  \item{actual}{
    adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
  }
  

  \item{cut   }{ number of cut points.}
  \item{smooth}{ logical: for smoothing the ROC curve.} 
  \item{calibrate}{If \code{TRUE}, compute the value of AUC by taking the level of the probabilities into account.}
  \item{linetype  }{ specification for the default plotting line type.} 
  \item{color}{ specification for the default plotting color.} 
  \item{size  }{ specification for the default plotting line size.} 
  \item{main  }{ overall title for the plot.} 
  \item{xlab  }{ title for the x axis.} 
  \item{ylab  }{ title for the y axis.} 
  \item{legend }{ logical: for adding legend to the ROC plot.} 
  \item{legend.size  }{ title for the x axis.} 
  \item{legend.position  }{ title for the y axis.} 
  \item{labels }{ for legends of the legend to the ROC plot.} 
  \item{auc}{ logical: to report AUC with legend.} 
  \item{theme}{ theme for the plot from the function \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}} of \code{R} package \code{\link[ggplot2:ggplot]{ggplot2::ggplot()}}.} 
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl}; Lucas Vogels \email{l.f.o.vogels@uva.nl} }

\seealso{ \code{\link{roc}}, \code{\link[pROC:plot.roc]{pROC::plot.roc()}}, \code{\link[pROC:auc]{pROC::auc()}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{compare}} }

\examples{
\dontrun{
# To generate multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 30, p = 6, size = 7, vis = TRUE )
   
# To Run sampling algorithm
bdgraph.obj <- bdgraph( data = data.sim, iter = 10000 )

# To compare the results
plotroc( bdgraph.ob2j, data.sim )
   
# To compare the results based on CGGMs approach
bdgraph.obj2 <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )

# To Compare the resultss
plotroc( list( bdgraph.obj, bdgraph.obj2 ), data.sim, legend = FALSE )
}
}

\keyword{hplot}

