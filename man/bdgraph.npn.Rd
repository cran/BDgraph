\name{bdgraph.npn}
\alias{bdgraph.npn}

\title{Nonparametric transfer}

\description{
This function transfer non-Gaussian data to Gaussian models.  
}
\usage{
bdgraph.npn(data, npn = "shrinkage", npn.thresh = NULL)
}

\arguments{
  \item{data}{ it could be a matrix or a data.frame of the data ( \eqn{n \times p}{n x p}). }
  \item{npn}{ 
  a character with three options "shrinkage" (default), "truncation", and "skeptic". Option "shrinkage" is for the 
  shrunken transformation, option "truncation" is for the truncated transformation, option "skeptic" is for the 
  non-paranormal skeptic transformation. For more details see references.
  }
  \item{npn.thresh}{ the truncation threshold; it is only for the truncated transformation (when
  \code{npn = "truncation"}). The default value is \eqn{1/(4n^{1/4} \sqrt{\pi \log(n)})}. }
}
\value{
\item{data}{ a \eqn{n \times p} matrix of transferred data, if \code{npn = "shrinkage" or "truncation"}, and 
  a non-paranormal correlation matrix (\eqn{p \times p}), if \code{npn = "skeptic"}.  
  }
}
\references{
Liu, H., Han, F., Yuan, M., Lafferty, J., and Wasserman, L., (2012). High Dimensional Semiparametric Gaussian Copula 
Graphical Models. \emph{Annals of Statistics} 

Zhao, T. and H. Liu. (2012). The huge Package for High-dimensional Undirected Graph Estimation in R. 
\emph{Journal of Machine Learning Research}

Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}

\author{ Abdolreza Mohammadi and Ernst Wit \cr }

\seealso{ \code{\link{bdgraph.sim}} and \code{\link{bdgraph}} }

\examples{
\dontrun{
  # generating the data from a random graph
  data.sim <- bdgraph.sim(n = 6, p = 4, size = 4)
  data     <- (data.sim $ data -3) ^ 4
  data
  
  # transfer the data by truncation 
  bdgraph.npn(data, npn = "truncation")
  
  # transfer the data by shrunken 
  bdgraph.npn(data, npn = "shrunken")

  # transfer the data by skeptic 
  bdgraph.npn(data, npn = "skeptic")
  }
}
