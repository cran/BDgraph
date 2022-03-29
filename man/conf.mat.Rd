\name{conf.mat}
\alias{conf.mat}

\title{	Confusion Matrix }
\description{ Create a Confusion Matrix. }

\usage{ 
conf.mat( pred, actual, cutoff = 0.5, proportion = FALSE, 
                     dnn = c( "Prediction", "Actual" ), ... )
}

\arguments{
    \item{pred  }{
        adjacency matrix corresponding to an estimated graph. 
        It can be an object with \code{S3} class "\code{bdgraph}" from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
    }
    
    \item{actual}{ adjacency matrix corresponding to the actual graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
	It can be a factor, numeric or character vector of responses (true class), typically encoded with 0 (controls) and 1 (cases).
    Only two classes can be used in a ROC curve. }
    \item{cutoff}{cutoff value for the case that \code{pred} is vector of probabilites. The default is 0.5. }  
	\item{proportion}{
		logical: FALSE (default) for a confusion matrix with number of cases. 
		         TRUE, for a confusion matrix with the proportion of cases.
	}
    \item{dnn   }{names to be given to the dimensions in the result (the dimnames names). }
    \item{\dots }{further arguments to be passed to \code{table}. }
}

\value{ the results of \code{table} on \code{pred} and \code{actual}. }

\author{Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{conf.mat.plot}}, \code{\link{compare}}, \code{\link{roc}}, \code{\link{bdgraph}} }

\examples{
\dontrun{
set.seed( 100 )

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

# Running sampling algorithm based on GGMs 
sample.ggm <- bdgraph( data = data.sim, method = "ggm", iter = 10000 )

# Confusion Matrix for GGM method
conf.mat( pred = sample.ggm, actual = data.sim )
}
}

\keyword{structure learning}

