\name{compare}
\alias{compare}

\title{ Graph structure comparison }

\description{ This function provides several measures to assess the performance of the graphical structure learning. }

\usage{ compare( pred, actual, main = NULL, vis = FALSE ) }

\arguments{
  \item{pred}{
        adjacency matrix corresponding to an estimated graph. 
        It can be an object with \code{S3} class "\code{bdgraph}" from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be an object of \code{S3} class \code{"select"}, from the function \code{\link[huge]{huge.select}} of \code{R} package \code{\link[huge]{huge}}.
        It also can be a \code{list} of above objects for comparing two or more different approaches.
        }

  \item{actual}{
    adjacency matrix corresponding to the true graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
	}
 
  \item{main}{ character vector giving the names for the result table. }
  \item{vis }{ logical: if TRUE, visualize the true graph and estimated graph structures. } 
}

\value{
	\item{True positive }{ number of correctly estimated links.}
	\item{True negative }{ number of true non-existing links which is correctly estimated.}
	\item{False positive}{ number of links which they are not in the true graph, but are incorrectly estimated.}
	\item{False negative}{ number of links which they are in the true graph, but are not estimated.}
	\item{F1-score      }{ weighted average of the \code{"positive predictive"} and \code{"true positive rate"}. The F1-score value reaches its best value at 1 and worst score at 0.}
	\item{Specificity   }{ Specificity value reaches its best value at 1 and worst score at 0.}
	\item{Sensitivity   }{ Sensitivity value reaches its best value at 1 and worst score at 0.}
	\item{MCC           }{ Matthews Correlation Coefficients (MCC) value reaches its best value at 1 and worst score at 0.}
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 
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
compare( sample.ggm, data.sim, main = c( "True", "GGM" ), vis = TRUE )
      
# Running sampling algorithm based on GCGMs
sample.gcgm <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )

# Comparing GGM and GCGM methods
compare( list( sample.ggm, sample.gcgm ), data.sim, 
         main = c( "True", "GGM", "GCGM" ), vis = TRUE )
}
}

\keyword{structure learning}

