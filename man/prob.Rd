\name{prob}
\alias{prob}

\title{
Posterior probabilities for most highest graphs or specific graph
}
\description{
According to the output of the BD-MCMC algorithm, this function gives us the 
posterior probability for the most highest graphs or a specific graph.
}
\usage{
prob(output, g = 2, A = NULL)
}

\arguments{
  \item{output}{
a list which is the result of BD-MCMC algorithm from the 'bdgraph' function. 
}
  \item{g}{
you can select the number of graphical models with highest probabilities (default is 1). 
}
  \item{A}{
upper triangular matrix in which \eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}. 
}
}
\value{
\item{list.A }{ a list which includes the grpahs with the highest posterior probabilities.}
\item{prob.A }{ a vector which includes posterior probabilities of all the graphs in 'list.A'.}
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
  data <- bdgraph.sim(n = 100, p = 8, vis = T)
  
  output <- bdgraph(data, meanzero = T, iter = 2000)
  
  prob(output, A = as.matrix(data $ A))
  
  prob(output, g = 3)
  }
}
