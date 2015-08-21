\name{plotroc}
\alias{plotroc}

\title{ ROC plot }

\description{
Draws the ROC curve according to the true graph structure for the object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}.
}

\usage{ plotroc( G, prob, prob2 = NULL, cut = 20, smooth = FALSE ) }

\arguments{
  \item{G}{
	The adjacency matrix corresponding to the true graph structure in which \eqn{g_{ij}=1} if there is a link between notes \eqn{i}{i} and
	\eqn{j}{j}, otherwise \eqn{g_{ij}=0}. It also can be an object of \code{S3} class \code{"sim"}, from function \code{\link{bdgraph.sim}}.
  }
  
  \item{prob}{ Upper triangular matrix which shows the stimated posterior probabilities for all possible links. It also can be an object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  
  \item{prob2}{ This option is for comparing two different estimations. Upper triangular matrix which shows the stimated posterior inclusion probabilities for all possible links. It also can be an object of \code{S3} class \code{"bdgraph"}, from function \code{\link{bdgraph}}. }
  
  \item{cut}{ Number of cut points. The default value is 20. }
  
  \item{smooth}{ Logical: for smoothing the ROC curve. The default is FALSE.} 
}

\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 
}

\author{Abdolreza Mohammadi and Ernst Wit \cr}

\seealso{ \code{\link{bdgraph}} and \code{\link{compare}} }

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

# runing sampling algorithm
output <- bdgraph( data = data.sim, iter = 10000 )
# comparing the result
plotroc( data.sim, output )

# To compare the result based on CGGMs approach
output2 <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )
# comparing the results
plotroc( data.sim, output, output2 )
legend( "bottomright", c("Gaussian", "Copula"), lty = c(1,2), col = c(1, 4))   
}
}
