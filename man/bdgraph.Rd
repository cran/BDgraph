\name{bdgraph}
\alias{bdgraph}

\title{
Graph selection based on birth-death MCMC algorithm 
}
\description{
This function is the BD-MCMC algorithm for model selection in undirected Gaussian 
graphical models based on birth-death MCMC methodology. It works for both low 
and high dimensional graphs.
}
\usage{
bdgraph(data, n = NULL, npn = "normal", mean = NULL, method = NULL, 
        g.prior = "Uniform", iter = 5000, b = 3, burnin = floor(iter / 2), 
	    thin = 1, lambda = NULL, D = NULL, g.start = "full", K.start = NULL, 
			         mc = 10, trace = TRUE, save.all = FALSE)
}

\arguments{
  \item{data}{
it could be a matrix or a data.frame of the data ( \eqn{n \times p}) or it could be a
covariance matrix as \eqn{S=x'x} which \eqn{x} is the data matrix.
}
  \item{n}{the number of observations. We need it if the "data" is a covariance matrix.}
  \item{npn}{ 
a character with four options "normal" (default), "shrinkage", "truncation", and "skeptic". Option "normal" means data 
are coming from multivariate normal distribution. Option "shrinkage" is for the shrunken transformation, option 
"truncation" is for the truncated transformation, option "skeptic" is for the non-paranormal skeptic transformation. 
For more details see references.
}
  \item{mean}{If \eqn{mean = 0} means the data have zero mean, otherwise the mean of the data should not be zero.}
  
  \item{method}{
a character with two options "fast" and "mc" to determine a type of BD-MCMC algorithm. Option "fast" is for a 
high-dimensional graphs (roughly graph with more than 8 nodes which is the default) and option "mc" is based on Monte Carlo
approximation for ratio of normalizing constant (for more details see reference paper).
}

  \item{g.prior}{
a character for selecting a prior distribution for the graph. It can be either "Uniform" (default) or "Poisson".
"Uniform" means discrete uniform distribution for prior distribution of the graph 
(\eqn{G \sim DU(\mathcal{G})} in which \eqn{\mathcal{G}} is all possible graphical models).
Also, "Poisson" means prior distribution for number of edges in the graph G has a truncated Poisson distribution
(\eqn{degree(G) \sim TP(\lambda)}.
}

  \item{iter}{the number of iteration for the BD-MCMC algorithm.}
  
  \item{b}{
the degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}. G-Wishart distribution is the prior
distribution of precision matrix.
}  
  
  \item{burnin}{the number of burn-in iteration for the BD-MCMC algorithm.}
  \item{thin}{option for regularly saving part of the iteration. Default is 1 (for saving all iteration).}

  \item{lambda}{a rate for prior distribution of graph when \eqn{degree(G) \sim TP(\lambda).}}

  \item{D}{
the positive definite matrix for G-Wishart distribution, \eqn{W_G(b,D)}. G-Wishart distribution is the prior 
distribution of precision matrix.
}
  \item{g.start}{
it shows the initial graph for the algorithm. It can be either "full" (default), "empty", "glasso", "mb", or "ct". 
"full" means the initial graph is a full graph and "empty" means a empty graph. "glasso" means 
the initial graph is a suggested graph by graphical lasso (glasso) method. "mb" means the initial graph is a 
suggested graph by Meinshausen-Buhlmann graph estimation (mb). "ct" means the initial graph is a 
suggested graph by correlation thresholding graph estimation (ct). We can use this option for high-dimensional cases 
which we need less iteration to converge, because we start with the graph which is almost near to the true graph. 
It also could be an object with S3 class "bdgraph"; with this option we can run the bdgraph from last objects of 
previous bdgraph run(see examples).     
} 

  \item{K.start}{it is a starting point for precision matrix. It has to be the positive definite matrix.}
 
  \item{mc}{
the number of iteration for Monte Carlo approximation of normalizing constant in G-Wishart distribution. 
We need it only when number of nodes in the graph is less than 8.
}

  \item{trace}{logical: if TRUE (default), prints out algorithm progress. }
  
  \item{save.all}{
logical: if FALSE (default), we save the adjacency matrices after burn-in. If TRUE, we save all matrices from 
starting point.}
}

\value{
an object with S3 class "bdgraph" as follow:
\item{sample.G}{a vector which includes the adjacency matrices for all iteration after burn-in.}

\item{weights}{a vector which includes the waiting times for all iteration after burn-in.}

\item{all.G}{a vector which includes the adjacency matrices for all iteration (includes burn-in iteration). 
We need it for checking the convergency of the BD-MCMC algorithm.}

\item{Khat}{estimation for precision matrix which is a mean of all samples from precision matrices.}
}

\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{ \code{\link{bdgraph.sim}} }

\examples{
\dontrun{
  # generating the data from a random graph
  data.sim <- bdgraph.sim(n = 100, p = 8, size = 10, vis = TRUE)
  
  # Comparing our result with true graph
  output <- bdgraph(data.sim, mean = 0, iter = 2000)
  
  compare(data.sim, output, colnames = c("True graph", "BDgraph"))
  
  output2 <- bdgraph(data.sim, mean = 0, iter = 2000, g.start = output)
  
  compare(data.sim, output, output2, colnames = c("True graph", "Frist run", "Second run"))
  }
}

