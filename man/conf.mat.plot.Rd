\name{conf.mat.plot}
\alias{conf.mat.plot}

\title{	Plot Confusion Matrix }
\description{ Plot a Confusion Matrix. }

\usage{ 
conf.mat.plot( pred, actual, cutoff = 0.5, conf.level = 0, margin = 1, 
                          color = c( "#ff83a8", "#83ff9b" ), ... )
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
    \item{conf.level}{confidence level used for the confidence rings on the odds ratios. Must be a single nonnegative number less than 1; if set to 0 (the default), confidence rings are suppressed. }  
    \item{margin}{numeric vector with the margins to equate. Must be one of 1 (the default), 2, or c(1, 2), which corresponds to standardizing the row, column, or both margins in each 2 by 2 table. Only used if std equals "margins". }  
    \item{color}{vector of length 2 specifying the colors to use for the smaller and larger diagonals of each 2 by 2 table. }  
    \item{\dots }{options to be passed to \code{fourfoldplot}. }
}

\author{Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link{conf.mat}}, \code{\link{compare}}, \code{\link{roc}}, \code{\link{bdgraph}} }

\examples{
\dontrun{
set.seed( 100 )

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

# Running sampling algorithm based on GGMs 
sample.ggm <- bdgraph( data = data.sim, method = "ggm", iter = 10000 )

# Confusion Matrix for GGM method
conf.mat.plot( pred = sample.ggm, actual = data.sim )
}
}

\keyword{structure learning}

