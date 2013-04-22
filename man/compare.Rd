\name{compare}
\alias{compare}

\title{
Comparing the results according to the true graph
}
\description{
With this function we can compare the result from BDMCMC algorithm or other approches, according to the true graph structure. 
}
\usage{
compare(true.g, est.g, est.g2 = NULL, colnames = NULL) 
}

\arguments{
  \item{true.g}{
adjacency for the true graph in which \eqn{true.g_{ij}=1}{true.gij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{true.g_{ij}=0}{true.gij = 0}. 
}
  \item{est.g}{
adjacency matrix for estimated graph in which \eqn{est.g_{ij}=1}{est.gij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{est.g_{ij}=0}{est.gij = 0}.}
  \item{est.g2}{
adjacency matrix for second estimated graph in which \eqn{est.g2_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{est.g2_{ij}=0}{est.g2ij = 0}.}
  \item{colnames}{
a character vector giving the column names for the result table.}
}
\value{
Returns an output of the BDMCMC algorithm which is list-like and contains the following:
\item{true positive}{the number of correctly estimated links.}
\item{true negative}{the number of true non-existing links which is correctly estimated.}
\item{false positive}{the number of links which they are not in the true graph, but are incorrectly estimated.}
\item{false negative}{the number of links which they are in the true graph, but are not estimated.}
\item{accuracy}{the number of true results (both true positives and true negatives) divided by the total number of
true and false results.}
\item{balanced F-score}{a weighted average of the "positive predictive value" and "true positive rate". F-score value 
reaches its best value at 1 and worst score at 0.}
\item{positive predictive value}{the number of correctly estimated links divided by the total number of links in the 
estimated graph.}
\item{true positive rate}{the number of correctly estimated links divided by the total number of links in the true 
graph.}
\item{false positive rate}{the false positive value divided by the total number of links in the true graph.}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{
Antonio Abbruzzo, Ivan Vujacic, Abdolreza Mohammadi, and Ernst Wit \cr
Maintainer: Abdolreza Mohammadi<a.mohammadi@rug.nl>
}
\seealso{
\code{\link{bdgraph}} and \code{\link{select}}
}
\examples{
\dontrun{
  # generate the data (100 observations) from multivariate normal distribution
  data <- bdgraph.sim(n = 100, p = 8, vis = T)
  
  # selecting the best graph according to BDMCMC algorithm
  output <- bdgraph(data $ data, meanzero = T, iter = 5000)
  bdgraph.g <- select(output)

  # selecting the best graph according to glasso by using huge package
  glasso <- huge(data $ data, method = "glasso")
  glasso <- huge.select(glasso, criterion = "stars")    
  
  # comparing the result by using "compare" function
  compare(true.g = data $ A, est.g = bdgraph.g, est.g2 = glasso $ refit, 
          colnames = c("true graph", "BDgraph", "glasso"))
  }
}
