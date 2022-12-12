\name{auc}
\alias{auc}

\title{ Compute the area under the ROC curve }

\description{ 
  This function computes the numeric value of area under the ROC curve (AUC) specifically for graph structure learning.
}

\usage{ 
auc( pred, actual, cut = 200, calibrate = TRUE )
}

\arguments{
    \item{pred  }{
        adjacency matrix corresponding to an estimated graph. 
        It can be an object with \code{S3} class "\code{bdgraph}" from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be a \code{\link{numeric}} or \code{\link{ordered}} vector of the same length than \code{actual}, containing the predicted
  value of each observation.
    }
    
    \item{actual}{ adjacency matrix corresponding to the actual graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class "\code{sim}" from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
	It can be a factor, numeric or character vector of responses (true class), typically encoded with 0 (controls) and 1 (cases).
    Only two classes can be used in a ROC curve. }
  \item{cut   }{ number of cut points.}
  \item{calibrate}{If \code{TRUE}, compute the value of AUC by taking the level of the probabilities into account.}
}

\value{
The numeric AUC value}

\references{
  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874, \doi{10.1016/j.patrec.2005.10.010}
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77, 
  \doi{10.1186/1471-2105-12-77}.
}

\author{Reza Mohammadi \email{a.mohammadi@uva.nl}; Lucas Vogels \email{l.f.o.vogels@uva.nl}}

\seealso{ \code{\link{plotroc}}, \code{\link[pROC:plot.roc]{pROC::plot.roc()}}, \code{\link[pROC:auc]{pROC::auc()}}, \code{\link[pROC:print.roc]{pROC::print.roc()}}, \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}}, \code{\link{compare}} }

\examples{
\dontrun{
set.seed( 5 )

# Generating multivariate normal data from a 'scale-free' graph
data.sim = bdgraph.sim( n = 200, p = 15, graph = "scale-free", vis = TRUE )

# Running BDMCMC algorithm 
sample.bdmcmc = bdgraph( data = data.sim, algorithm = "bdmcmc", iter = 10000 )

BDgraph::auc( pred = sample.bdmcmc, actual = data.sim, calibrate = TRUE )

}
}

\keyword{structure learning}

