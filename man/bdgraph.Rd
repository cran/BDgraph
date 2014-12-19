\name{bdgraph}
\alias{bdgraph}

\title{
	Birth-death MCMC algorithm for graphical models 
}
\description{
	This function is the BD-MCMC algorithm for model selection in undirected Gaussian 
	graphical models based on birth-death MCMC methodology. 
}
\usage{
bdgraph( data, n = NULL, method = "ggm", iter = 5000, burnin = iter / 2, 
         b = 3, D = NULL, Gstart = "empty" )
}

\arguments{
	\item{data}{
	It could be a matrix or a data.frame of the data ( \eqn{n \times p}) or a covariance matrix as \eqn{S=x'x} which \eqn{x} is the data matrix. It also could be an object of class "simulate".
	}

	\item{n}{The number of observations. It is needed if the "data" is a covariance matrix.}
	  
	\item{method}{
		A character with two options "ggm" (defult) and "cggm". 
		Option "ggm" is based on Gaussian graphical models which is based on gaussianty assumiton.
		Option "cggm" is based on copula Gaussian graphical models which is for the data that not follow gaussianty assumiton e.g. continuous non-Gaussian, discrete, or mixed data sets.
	}

	\item{iter}{The number of iteration for the BD-MCMC algorithm.}
	\item{burnin}{The number of burn-in iteration for the BD-MCMC algorithm.}

	\item{b}{
		The degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}. G-Wishart distribution is the prior
		distribution of precision matrix. The default is 3.
	}
	\item{D}{
		The positive definite matrix for G-Wishart distribution, \eqn{W_G(b,D)}. G-Wishart distribution is the prior 
		distribution of precision matrix. The default is an identity matrix.
	}
	\item{Gstart}{
		It is corresponding to a starting point for graph. It can be either "full" (default), "empty", or an object with S3 class "bdgraph". 
		Option "full" means the initial graph is a full graph and "empty" means a empty graph. 
		Gstart also could be an object with S3 class "bdgraph"; With this option we could run the bdgraph algorithm from last objects of 
		previous bdgraph run (see examples).     
	} 
}

\value{
	an object with S3 class "bdgraph" as follow:

	\item{sampleGraphs}{A vector which includes the adjacency matrices for all iteration after burn-in.}
	\item{graphWeights}{A vector which includes the waiting times for all iteration after burn-in.}

	\item{allGraphs}{A vector which includes the adjacency matrices for all iteration (includes burn-in iteration). 
		  It is needed to check the convergency of the BD-MCMC algorithm.}

	\item{allWeights}{A vector which includes the waiting times for all iteration (includes burn-in iteration).}

	\item{Khat}{Estimation for precision matrix which is a mean of all samples from precision matrices.}
}

\references{
	Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
	graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}

\author{ Abdolreza Mohammadi and Ernst Wit }

\seealso{ \code{\link{bdgraph.sim}} }

\examples{
\dontrun{
# generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 20, p = 6, size = 7, vis = TRUE )

output <- bdgraph( data = data.sim, iter = 1000 )

summary(output)

# To compare our result with true graph
compare( data.sim, output, colnames = c("True graph", "BDgraph") )

output2 <- bdgraph( data = data.sim, iter = 5000, Gstart = output )

compare( data.sim, output, output2, colnames = c("True graph", "Frist run", "Second run") )

# generating mixed data from a 'scale-free' graph
data.sim <- bdgraph.sim( n = 100, p = 6, type = "mixed", graph = "scale-free", vis = TRUE )

output <- bdgraph( data = data.sim, method = "cggm", iter = 10000 )

summary( output )

compare( data.sim, output )	  
}
}

