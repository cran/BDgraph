\name{compare}
\alias{compare}

\title{ Comparing the result }

\description{ With this function, we can check the performance of our method and compare it with other alternative approaches. }

\usage{ compare( G, est, est2 = NULL, est3 = NULL, colnames = NULL, vis = FALSE ) }

\arguments{
  \item{G}{
	The adjacency matrix corresponding to the true graph in which \eqn{g_{ij}=1} if there is a link between notes \eqn{i}{i} and
	\eqn{j}{j}, otherwise \eqn{g_{ij}=0}. It also can be an object with \code{S3} class \code{"simulate"} from function \code{\link{bdgraph.sim}}.
	}
  \item{est}{ An adjacency matrix corresponding to an estimated graph. It also can be an object with \code{S3} class \code{"bdgraph"} from function \code{\link{bdgraph}} or \code{"select"} (from \pkg{huge} package). }
  \item{est2}{ An adjacency matrix corresponding to an estimated graph. It also can be an object with \code{S3} class \code{"bdgraph"} from function \code{\link{bdgraph}} or "select" (from \pkg{huge} package). }
  \item{est3}{ An adjacency matrix corresponding to an estimated graph. It also can be an object with \code{S3} class \code{"bdgraph"} from function \code{\link{bdgraph}} or "select" (from \pkg{huge} package). }
  
  \item{colnames}{ A character vector giving the column names for the result table. }
  
  \item{vis}{ Logical: if TRUE you will see a plot result. The default is FALSE. } 
}

\value{
	\item{True positive}{ The number of correctly estimated links.}
	\item{True negative}{ The number of true non-existing links which is correctly estimated.}
	\item{False positive}{ The number of links which they are not in the true graph, but are incorrectly estimated.}
	\item{False negative}{ The number of links which they are in the true graph, but are not estimated.}
	\item{Accuracy}{ the number of true results (both true positives and true negatives) divided by the total number of true and false results.}
	\item{Balanced F-score}{ A weighted average of the \code{"positive predictive"} and \code{"true positive rate"}. F-score value reaches its best value at 1 and worst score at 0.}
	\item{Positive predictive}{ The number of correctly estimated links divided by the total number of links in the estimated graph.}
	\item{True positive rate}{ The number of correctly estimated links divided by the total number of links in the true graph.}
	\item{False positive rate}{ The false positive value divided by the total number of links in the true graph.}
}
\references{
Mohammadi, A. and E. Wit (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. and E. Wit (2015). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Arxiv preprint arXiv:1501.05108v2} 

Mohammadi, A., F. Abegaz Yazew, E. van den Heuvel, and E. Wit (2015). Bayesian Modeling of Dupuytren Disease Using Gaussian Copula Graphical Models, \emph{Arxiv preprint arXiv:1501.04849v2} 
}
\author{Abdolreza Mohammadi, Antonio Abbruzzo, Ivan Vujacic, and Ernst Wit \cr}

\seealso{ \code{\link{bdgraph}} and \code{\link{select}} }

\examples{
\dontrun{
# generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

# running sampling algorithm
output <- bdgraph( data = data.sim, iter = 10000 )
# comparing the result
compare( data.sim, output, colnames = c("true", "BDgraph"), vis = TRUE )

# estimate the true graph by 'huge' package
require( huge )  
huge.g <- huge( data.sim $ data, method = "mb" )
huge.g <- huge.select( huge.g )    

# comparing the result
compare( data.sim, output, huge.g, colnames = c("true", "BDgraph", "huge") , vis = TRUE )
}
}
