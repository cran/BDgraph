\name{compare}
\alias{compare}

\title{
Comparing the results according to the true graph
}
\description{
With this function we can compare the result from BD-MCMC algorithm or other approches, according to the true graph structure. 
}
\usage{
compare(true, estimate, estimate2 = NULL, colnames = NULL, plot = TRUE) 
}

\arguments{
  \item{true}{
adjacency matrix for the true graph in which \eqn{true_{ij}=1}{trueij = 1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{true_{ij}=0}{trueij = 0}. 
}
  \item{estimate}{ adjacency matrix for estimated graph. }
  
  \item{estimate2}{ adjacency matrix for second estimated graph. }
  
  \item{colnames}{ a character vector giving the column names for the result table. }
  
  \item{plot}{ logical: if TRUE (default) you will see the plot result.} 
  }
\value{
\item{true positive}{the number of correctly estimated links.}
\item{true negative}{the number of true non-existing links which is correctly estimated.}
\item{false positive}{the number of links which they are not in the true graph, but are incorrectly estimated.}
\item{false negative}{the number of links which they are in the true graph, but are not estimated.}
\item{accuracy}{the number of true results (both true positives and true negatives) divided by the total number of
true and false results.}
\item{balanced F-score}{a weighted average of the "positive predictive" and "true positive rate". F-score value 
reaches its best value at 1 and worst score at 0.}
\item{positive predictive}{the number of correctly estimated links divided by the total number of links in the 
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
}
\seealso{
\code{\link{bdgraph}} and \code{\link{select}}
}
\examples{
\dontrun{
  # generate the data from multivariate normal distribution
  data.sim <- bdgraph.sim(n = 100, p = 8, size = 9, vis = T)
  
  # selecting the best graph according to BD-MCMC algorithm
  output <- bdgraph(data.sim, meanzero = T, iter = 5000)
  compare(data.sim, output, colnames = c("true graph", "BDgraph"))

  # selecting the best graph by huge package
  huge.g <- huge(data.sim $ data, method = "mb")
  huge.g <- huge.select(huge.g)    
  
  # comparing the result by using "compare" function
  compare(data.sim, output, estimate2 = huge.g $ refit, colnames = c("true graph", "BDgraph", "huge"))
  }
}
