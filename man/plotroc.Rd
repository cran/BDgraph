\name{plotroc}
\alias{plotroc}

\title{ ROC plot }

\description{
Plot the ROC curve according to the true graph structure for the "bdgraph" output.
}

\usage{
plotroc( G, prob, prob2 = NULL, cut = 20, smooth = FALSE )
}

\arguments{
  \item{G}{
The adjacency matrix corresponding to the true graph structure in which \eqn{G_{ij}=1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{G_{ij}=0}. It also can be an object with S3 class "simulate".
}
  \item{prob}{ Upper triangular matrix which shows the stimated posterior inclusion probabilities for all possible edges. It also can be an object with S3 class "bdgraph". }
  
  \item{prob2}{ This option is for comparing two different estimations. Upper triangular matrix which shows the stimated posterior inclusion probabilities for all possible edges. It also can be an object with S3 class "bdgraph". }
  
  \item{cut}{ Number of cut point. Default is 20. }
  
  \item{smooth}{ Logical: for smoothing the ROC curve. The default is FALSE.} 
  }

\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit \cr}

\seealso{ \code{\link{bdgraph}} and \code{\link{compare}} }

\examples{
\dontrun{
  # generating synthetic multivariate normal data from a 'random' graph
  data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
  
  # selecting the best graph according to BD-MCMC algorithm
  output <- bdgraph( data = data.sim, iter = 10000 )
  # comparing the result
  plotroc( data.sim, output )
  
  plotroc( data.sim, output, smooth = TRUE )
  }
}
