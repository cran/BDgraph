\name{bdgraph.sim}
\alias{bdgraph.sim}

\title{
Synthetic graph data generator
}
\description{
Implements a synthetic graph data generation with different graph structures, including 
"random", "cluster", "scale-free", "hub" "fixed", and "circle". The default is "random".
Based on the graph structure, it generate four differents data sets which is
multivariate Gaussian, non-Gaussian, discrete, or mixed data.  
}
\usage{
bdgraph.sim( n = 2, p = 10, graph = "random", size = NULL, prob = 0.2, 
             class = NULL, type = "Gaussian", cut = 4, b = 3, D = diag(p), 
             K = NULL, sigma = NULL, mean = 0, vis = FALSE )
}

\arguments{
  \item{n}{The number of samples required. The default value is 2.}
  \item{p}{The number of variables (nodes). The default value is 10.}
  \item{graph}{The graph structure with option "random" (as a default), "cluster", "scale-free", "hub" "fixed", and "circle". 
     It also could be adjacency matrix corresponding to the graph structure (an upper triangular matrix in which 
     \eqn{g_{ij}=1} if there is a link between notes \eqn{i} and \eqn{j}, otherwise \eqn{g_{ij}=0}). }
  \item{size}{The number of links in true graph (graph size).}
  \item{prob}{For "random" graph, it is the probability that a pair of nodes has a link. The default value is 0.2.}
  \item{class}{Is number of the classes for "cluster" option.}
  
  \item{type}{Type of data which has four options "Gaussian" (as a default), "non-Gaussian", "discrete", and "mixed".
  For option "Gaussian", data is generated from multivariate normal distribution.
  For option "non-Gaussian", data is transfered multivariate normal distribution to continuous multivariate non-Gaussian distribution.
  For option "discrete", data is transfered from multivariate normal distribution to discrete multivariat distribution.
  For option "mixed", data is transfered from multivariate normal distribution to mixture of 'count', 'ordinal', 'non-Gaussian', 'binary' and 'Gaussian', respectively. 
  }
  
  \item{cut}{ For simulated discrete data (type = "discrete"), it is number of categories for discrete variables. The default value is 4.}

  \item{b}{The degree of freedom for G-Wishart distribution, \eqn{W_G(b, D)}. The default is 3.}
  \item{D}{The positive definite \eqn{(p \times p)} "scale" matrix for G-Wishart distribution, \eqn{W_G(b, D)}. The default is an identity matrix.}

  \item{K}{A positive-definite symmetric matrix specifying the precision matrix of the variables. It is for "fixed" 
  option of graph. }
  \item{sigma}{A positive-definite symmetric matrix specifying the covariance matrix of the variables. It is for 
  "fixed" option of graph.}
  
  \item{mean}{A vector specifying the mean of the variables. The default value is zero.}
  \item{vis}{visualize the true graph pattern. The default value is FALSE.}
}
\value{
an object with S3 class "simulate" as follow:
\item{data}{Generated data which is a \eqn{n \times p}{n x p} matrix.}
\item{sigma}{The covariance matrix of the generated data.}
\item{K}{The precision matrix of the generated data.}
\item{G}{The adjacency matrix corresponding to the true graph structure.}
}
\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi and Ernst Wit}
\seealso{\code{\link{bdgraph}}}

\examples{
	\dontrun{
		# generating synthetic multivariate normal data from a 'random' graph
		data.sim <- bdgraph.sim( n = 50, p = 10, prob = 0.3, vis = TRUE )
		print(data.sim)

		# generating synthetic multivariate normal data from a 'hub' graph
		data.sim <- bdgraph.sim( n = 3, p = 6, graph = "hub", vis = FALSE )
		data.sim $ data

		# generating synthetic mixed data from a 'hub' graph 
		data.sim <- bdgraph.sim( n = 10, p = 8, graph = "hub", type = "mixed" )
		data.sim $ data
	}
}
