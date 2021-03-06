\name{compare}
\alias{compare}

\title{ Graph structure comparison }

\description{ This function provides several measures to assess the performance of the graphical structure learning. }

\usage{ 
compare( target, est, est2 = NULL, est3 = NULL, est4 = NULL, main = NULL, 
         vis = FALSE ) 
}

\arguments{
  \item{target}{
    An adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class \code{"sim"} from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
	}
  \item{est,
        est2,
        est3,
        est4}{
        An adjacency matrix corresponding to an estimated graph. 
        It can be an object with \code{S3} class \code{"bdgraph"} from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be an object of \code{S3} class \code{"select"}, from the function \code{\link[huge]{huge.select}} of \code{R} package \code{\link[huge]{huge}}.
        Options \code{est2, est3} and \code{est4} are for comparing two or more different approaches.
        }
 
  \item{main}{ A character vector giving the names for the result table. }
  \item{vis }{ Visualize the true graph and estimated graph structures. } 
}

\value{
	\item{True positive }{ The number of correctly estimated links.}
	\item{True negative }{ The number of true non-existing links which is correctly estimated.}
	\item{False positive}{ The number of links which they are not in the true graph, but are incorrectly estimated.}
	\item{False negative}{ The number of links which they are in the true graph, but are not estimated.}
	\item{F1-score      }{ A weighted average of the \code{"positive predictive"} and \code{"true positive rate"}. The F1-score value reaches its best value at 1 and worst score at 0.}
	\item{Specificity   }{ The Specificity value reaches its best value at 1 and worst score at 0.}
	\item{Sensitivity   }{ The Sensitivity value reaches its best value at 1 and worst score at 0.}
	\item{MCC           }{ The Matthews Correlation Coefficients (MCC) value reaches its best value at 1 and worst score at 0.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645 

Letac, G., Massam, H. and Mohammadi, R. (2018). The Ratio of Normalizing Constants for Bayesian Graphical Gaussian Model Selection, \emph{arXiv preprint arXiv:1706.04416v2} 

Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845
}
\author{Reza Mohammadi \email{a.mohammadi@uva.nl}, Antonio Abbruzzo, and Ivan Vujacic}

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{bdgraph.sim}}, \code{\link{plotroc}} }

\examples{
\dontrun{
# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )
    
# Running sampling algorithm based on GGMs 
sample.ggm <- bdgraph( data = data.sim, method = "ggm", iter = 10000 )
   
# Comparing the results
compare( data.sim, sample.ggm, main = c( "True", "GGM" ), vis = TRUE )
      
# Running sampling algorithm based on GCGMs
sample.gcgm <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )

# Comparing GGM and GCGM methods
compare( data.sim, sample.ggm, sample.gcgm, main = c( "True", "GGM", "GCGM" ), vis = TRUE )
}
}

\keyword{structure learning}

