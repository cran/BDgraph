\name{bdgraph.sim}
\alias{bdgraph.sim}

\title{
Data generator according to graph structure
}
\description{
Implements the data generation from multivariate normal distributions with different graph structures, including 
"random", "fixed". The default is "random".
}
\usage{
bdgraph.sim(n = 1, p = 10, graph = "random", prob = NULL, v = NULL, u = NULL, 
             A = NULL, K = NULL, sigma = NULL, vis = FALSE)
}

\arguments{
  \item{n}{
the number of samples required.
}
  \item{p}{
the number of variables (nodes). 
}
  \item{graph}{
the graph structure with option "random" and "fixed". The default is "random".
}
  \item{prob}{
for "random" graph, it is the probability that a pair of nodes has a link. The default value is 0.2.
}
  \item{v}{
the off-diagonal elements of the precision matrix, controlling the magnitude of partial correlations with u. 
The default value is 0.3.
}
  \item{u}{
A positive number being added to the diagonal elements of the precision matrix, to control the magnitude of partial 
correlations. The default value is 0.1.
}
  \item{A}{
upper triangular matrix in which \eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}.
}
  \item{K}{
a positive-de?nite symmetric matrix specifying the precision matrix of the variables. It is for "fixed" option of graph.
}
  \item{sigma}{
a positive-de?nite symmetric matrix specifying the covariance matrix of the variables. It is for "fixed" option of graph.
}
  \item{vis}{
visualize the graph pattern. The default value is FALSE.
}
}
\value{
\item{data}{the \eqn{n \times p}{n x p} matrix for the generated data}
\item{sigma}{the covariance matrix of the generated data}
\item{K}{the precision matrix of the generated data}
\item{A}{the adjacency matrix of true graph structure}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdgraph}}
}
\examples{
\dontrun{
  # generate the data (100 observations) from multivariate normal distribution
  data.sim <- bdgraph.sim(n = 100, p = 8, vis = T)
  data.sim $ data
  }
}