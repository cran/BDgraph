\name{compare}
\alias{compare}

\title{
Comparing the result
}
\description{
With this function we can compare the result from BD-MCMC algorithm or other approaches, usually based to the true graph structure. 
}
\usage{
compare( G, est, est2 = NULL, colnames = NULL, vis = FALSE ) 
}

\arguments{
  \item{G}{
The adjacency matrix corresponding to the true graph in which \eqn{G_{ij}=1} if there is a link between notes \eqn{i}{i} and
\eqn{j}{j}, otherwise \eqn{G_{ij}=0}. It also can be an object with S3 class "simulate".
}
  \item{est}{ Adjacency matrix corresponding to the estimated graph. It also can be an object with S3 class "bdgraph" or "select" (from huge package). }
  
  \item{est2}{ Adjacency matrix corresponding to the estimated graph. It also can be an object with S3 class "bdgraph" or "select" (from huge package). }
  
  \item{colnames}{ A character vector giving the column names for the result table. }
  
  \item{vis}{ Logical: if TRUE you will see the plot result. The default is FALSE.} 
  }
\value{
\item{True positive}{the number of correctly estimated links.}
\item{True negative}{the number of true non-existing links which is correctly estimated.}
\item{False positive}{the number of links which they are not in the true graph, but are incorrectly estimated.}
\item{False negative}{the number of links which they are in the true graph, but are not estimated.}
\item{Accuracy}{the number of true results (both true positives and true negatives) divided by the total number of
true and false results.}
\item{Balanced F-score}{a weighted average of the "positive predictive" and "true positive rate". F-score value 
reaches its best value at 1 and worst score at 0.}
\item{Positive predictive}{the number of correctly estimated links divided by the total number of links in the 
estimated graph.}
\item{True positive rate}{the number of correctly estimated links divided by the total number of links in the true 
graph.}
\item{False positive rate}{the false positive value divided by the total number of links in the true graph.}
}
\references{
Mohammadi, A. and Wit, E. C. (2014). Bayesian structure learning in sparse Gaussian 
graphical models, Bayesian Analysis, acceped. \url{http://arxiv.org/abs/1210.5371v6}
}
\author{Abdolreza Mohammadi, Antonio Abbruzzo, Ivan Vujacic, and Ernst Wit \cr}

\seealso{ \code{\link{bdgraph}} and \code{\link{select}} }

\examples{
\dontrun{
# generating synthetic multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

# selecting the best graph according to BD-MCMC algorithm
output <- bdgraph( data = data.sim, iter = 10000 )
# comparing the result
compare( data.sim, output, colnames = c("true", "BDgraph"), vis = TRUE )

# selecting the best graph by 'huge' package
require( huge )  
huge.g <- huge( data.sim $ data, method = "mb" )
huge.g <- huge.select( huge.g )    

# comparing the result
compare( data.sim, output, huge.g, colnames = c("true", "BDgraph", "huge") , vis = T )
}
}
