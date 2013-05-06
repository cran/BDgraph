\name{bdgraph}
\alias{bdgraph}

\title{
Graph determination based on BD-MCMC algorithm 
}
\description{
This function is the BD-MCMC algorithm for model selection in undirected Gaussian 
graphical models based on birth-death MCMC methodology. It works for both low 
and high dimensional graphs.
}
\usage{
bdgraph(data, n = NULL, meanzero = FALSE, model = NULL, iter = 5000, 
       burnin = floor(iter / 2), skip = 1, birthrate = NULL, prior.g = "Uniform", 
	       b = 3, D = NULL, start.g = "full", MCiter = 10, summary = FALSE, 
						       trace = TRUE, save.all = FALSE, last.objects = NULL)
}

\arguments{
  \item{data}{
it could be a matrix or a data.frame of the data ( \eqn{n \times p}{n x p}) or it could be a
covariance matrix as \eqn{S=x'x}{S = x' x} which \eqn{x}{x} is the data.
}
  \item{n}{
the number of observations.
}
  \item{meanzero}{
logical: if TRUE that means the data have zero mean; if FALSE (default) the mean of the data should not be zero.
}
  \item{model}{
a character with two options "high" and "low" to determine a type of BD-MCMC algorithm. Option "high" is for a 
high-dimensional graphs (roughly graph with more than 8 nodes) and option "low" is for a low-dimensional graphs.
For more ditails see reference paper.
}
  \item{iter}{
the number of iterations for the BD-MCMC algorithm.
}
  \item{burnin}{
the number of burn-in iterations for the BD-MCMC algorithm.
}
  \item{skip}{
it is the option for regularly saving part of the iterations. Default is 1 (for saving all iterations).
}
  \item{birthrate}{
the birth rates for the birth-death process.
}
  \item{prior.g}{
a character for selecting a prior distribution for the graph. It can be either "Uniform" (default) or "Poisson".
"Uniform" means discrete uniform distribution for prior distribution of the graph 
(\eqn{G \sim DU(\mathcal{G})}{G ~ DU(G)} in which \eqn{\mathcal{G}}{G} is all possible graphical models).
Also, "Poisson" means prior distribution for number of edges in the graph G has a truncated Poisson distribution
(\eqn{degree(G) \sim TP(\gamma_b)}{degree(G) ~ TP(birthrate)} in which for our model \eqn{\gamma_b}{gammab} equal 
with a 'birthrate' value).
}
  \item{b}{
the degree of freedom for G-Wishart distribution, \eqn{W_G(b,D)}{WG(b,D)}. G-Wishart distribution is the prior
distribution of precision matrix.
}
  \item{D}{
the positive definite matrix for G-Wishart distribution, \eqn{W_G(b,D)}{WG(b,D)}. G-Wishart distribution is the prior 
distribution of precision matrix.
}
  \item{start.g}{
it shows the starting graph for BD-MCMC algorithm. It is an upper triangular matrix in which 
\eqn{g_{ij}=1}{gij = 1} if there is a link between nodes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{g_{ij}=0}{gij = 0}.
It can be either "full" (default) or "empty" or "glasso". "full" means we start with full graph. "empty" means we start
with empty graph. "glasso" means we start with a graph which is suggested by huge package. We can use this option for 
high-dimentional cases which we need less iterations to converge, because we start with the graph which is almost near 
to the true graph.  
}  
  \item{MCiter}{
the number of iterations for Monte Carlo approximation of normalizing constant in G-Wishart distribution. 
We need it only when number of nodes in the graph is less than 8.
}
  \item{summary}{
logical: if TRUE you will see the summary result of BD-MCMC algorithm; if FALSE (default), the output will be a 
list as below.
}
  \item{trace}{
logical: if trace = TRUE, tracing number of iterations in BD-MCMC algorithm printing is abled. The default value is TRUE.
}
  \item{save.all}{
logical: if save.all = FALSE, we save the matrix A after burn-in iteration. The default value is FALSE. If save.all = FALSE,
we save all matrix A from starting point of iteration.
}
  \item{last.objects}{
If you would like to run the BD-MCMC algorithm for several times, with this obtion you can run the bdgraph from last 
object of previous run. 
}
}
\value{
an object with S3 class "bdgraph" as follow:
\item{sample.G}{a vector which includes the adjacency matrices for all iterations after burn-in.}
\item{lambda}{a vector which includes the waiting times for all iterations after burn-in.}
\item{all.G}{a vector which includes the adjacency matrices for all iterations (includes burn-in iterations). 
We need it for checking the convergency of the BD-MCMC algorithm.}
\item{Khat}{estimation for precision matrix which is a mean of all samples from precision matrices.}
\item{p}{the number of nodes.}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdgraph.sim}}
}
\examples{
\dontrun{
  # generate the data from multivariate normal distribution
  data.sim <- bdgraph.sim(n = 100, p = 8, size = 10, vis = T)
  
  bdgraph(data.sim, meanzero = T, iter = 2000, summary = T)
  
  # for compare our result with true graph, we can run
  output <- bdgraph(data.sim, meanzero = T, iter = 2000)
  
  compare(data.sim, output, colnames = c("True graph", "BDgraph"))
  
  # we also can run which talk more time
  bdgraph(data.sim, meanzero = T, model = "low", iter = 2000, summary = T)
  }
}

