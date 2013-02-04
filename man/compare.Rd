\name{compare}
\alias{compare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparing the result according to the true graph
}
\description{
According to the true graph structure, with "compare" function we can compare the result for bdmcmc algorithm. 
}
\usage{
compare(est.g, true.g) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est.g}{
adjacency matrix for estimated graph in which \eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}.}
  \item{true.g}{
adjacency for the true graph in which \eqn{a_{ij}=1}{aij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{a_{ij}=0}{aij = 0}. 
}
}
\value{
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
Antonio Abbruzzo, Ivan Vujacic, Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdmcmc}} and \code{\link{select.g}}
}
\examples{
\dontrun{
  require(MASS)
  
  p <- 8 # number of nodes 
  # "truK" is the precision matrix of true graph
  truK <- diag(p)
  for (i in 1:(p-1)) truK[i,i+1] <- truK[i+1,i] <- 0.5
  truK[1,p] <- truK[p,1] <- 0.4
  truK # precision matrix of the true graph
  
  # generate the data (200 observations) from multivariate normal 
  # distribution with mean zero and percision matrix "truK"
  data <- mvrnorm(200, c(rep(0,p)), solve(truK))  
  
  # selecting the best graph according to bdmcmc algorithm
  output <- bdmcmc(data, meanzero = T, iter = 5000)
  bdmcmc.g <- select.g(output)

  # selecting the best graph according to glasso by using huge package
  huge.g <- huge(data, method = "glasso")
  huge.g <- huge.select(huge.g, criterion = "stars")    
  
  true.g <- ceiling (truK)
  true.g # matrix A shows circle graph with 8 links and 8 nodes which is the true graphical model
  
  # comparing the result by using "compare" function
  compare.true <- compare(est.g = true.g, true.g = true.g)
  compare.bdmcmc <- compare(est.g = bdmcmc.g, true.g = true.g)
  compare.huge <- compare(est.g = huge.g$refit, true.g = true.g)  
  compare.all <- cbind(compare.true, compare.bdmcmc, compare.huge)
  colnames(compare.all) <- c("True graph", "BDgraph", "huge")
  compare.all   
  }
}
