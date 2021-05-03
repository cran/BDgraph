\name{roc}
\alias{roc}

\title{ Build a ROC curve }

\description{ 
  This function builds a ROC curve specifically for graph structure learning and returns a \dQuote{roc} object, a list of class
  \dQuote{roc}. This object can be \code{prin}ted, \code{plot}ted, or
  passed to the functions \code{\link[pROC:roc]{pROC::roc()}}, \code{\link[pROC:ci]{pROC::ci()}}, \code{\link[pROC:smooth.roc]{pROC::smooth.roc()}}
  and \code{\link[pROC:coords]{pROC::coords()}}. Additionally, two \code{roc} objects can be compared with \code{\link[pROC:roc.test]{pROC::roc.test()}}. 
  This function is based on the \code{\link[pROC]{roc}} function of \code{R} package \code{pROC}.
}

\usage{ 
roc( pred, actual, auc = TRUE, smooth = FALSE, plot = FALSE, ... )
}

\arguments{
    \item{pred  }{
        An adjacency matrix corresponding to an estimated graph. 
        It can be an object with \code{S3} class \code{"bdgraph"} from function \code{\link{bdgraph}}. 
        It can be an object of \code{S3} class \code{"ssgraph"}, from the function \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}} of \code{R} package \code{\link[ssgraph:ssgraph]{ssgraph::ssgraph()}}.
        It can be a \code{\link{numeric}} or \code{\link{ordered}} vector of the same length than \code{actual}, containing the predicted
  value of each observation.
    }
    
    \item{actual}{ An adjacency matrix corresponding to the actual graph structure in which \eqn{a_{ij}=1} if there is a link between notes \eqn{i}{i} and \eqn{j}{j}, otherwise \eqn{a_{ij}=0}.
	It can be an object with \code{S3} class \code{"sim"} from function \code{\link{bdgraph.sim}}.
	It can be an object with \code{S3} class \code{"graph"} from function \code{\link{graph.sim}}.
	It can be a factor, numeric or character vector of responses (true class), typically encoded with 0 (controls) and 1 (cases).
    Only two classes can be used in a ROC curve. }
 \item{smooth}{if TRUE, the ROC curve is passed to \code{\link{smooth}}
    to be smoothed.
  }
  \item{auc}{compute the area under the curve (AUC)? If \code{TRUE}
    (default), additional arguments can be passed to \code{\link[pROC:auc]{pROC::auc()}}.
  }
  
  \item{plot}{plot the ROC curve? If \code{TRUE}, additional
    arguments can be passed to \code{\link[pROC:plot.roc]{pROC::plot.roc()}}.
  }
  \item{\dots}{ further arguments to be passed to \code{\link[pROC:roc]{pROC::roc()}}. }
}

\value{
  If the data contained any \code{NA} value and \code{na.rm=FALSE}, \code{NA} is
  returned. Otherwise, if \code{smooth=FALSE}, a list of class
  \dQuote{roc} with the following fields: 
  \item{auc}{if called with \code{auc=TRUE}, a numeric of class \dQuote{auc} as
    defined in \code{\link[pROC:auc]{pROC::auc()}}.
  }
  \item{ci}{if called with \code{ci=TRUE}, a numeric of class \dQuote{ci} as
    defined in \code{\link[pROC:ci]{pROC::ci()}}.
  }
  \item{response}{the response vector. Patients whose response is not
    \code{\link{\%in\%}} \code{levels} are discarded. If \code{NA} values
    were removed, a \code{na.action} attribute similar
    to \code{\link{na.omit}} stores the row numbers. 
  }
  \item{predictor}{the predictor vector converted to numeric as used to build the ROC
    curve. Patients whose response is not \code{\link{\%in\%}} \code{levels} are discarded. If
    \code{NA} values were removed, a \code{na.action} attribute similar
    to \code{\link{na.omit}} stores the row numbers.
  }
  \item{original.predictor, original.response}{the response and predictor vectors as passed in argument.}
  \item{levels}{the levels of the response as defined in argument.}
  \item{controls}{the predictor values for the control observations.}
  \item{cases}{the predictor values for the cases.}
  \item{percent}{if the sensitivities, specificities and AUC are
    reported in percent, as defined in argument.
  }
  \item{direction}{the direction of the comparison, as defined in argument.}
  \item{fun.sesp}{the function used to compute sensitivities and specificities.
    Will be re-used in bootstrap operations.}
  \item{sensitivities}{the sensitivities defining the ROC curve.}
  \item{specificities}{the specificities defining the ROC curve.}
  \item{thresholds}{the thresholds at which the sensitivities and
    specificities were computed. See below for details.
  }
  \item{call}{how the function was called. See \code{\link{match.call}} for
    more details.
  }

  If \code{smooth=TRUE} a list of class \dQuote{smooth.roc} as returned
  by \code{\link[pROC:smooth]{pROC::smooth()}}, with or without additional elements
  \code{auc} and \code{ci} (according to the call).
}

\references{
  Tom Fawcett (2006) ``An introduction to ROC analysis''. \emph{Pattern
    Recognition Letters} \bold{27}, 861--874. DOI:
  \doi{10.1016/j.patrec.2005.10.010}.
  
  Xavier Robin, Natacha Turck, Alexandre Hainard, \emph{et al.}
  (2011) ``pROC: an open-source package for R and S+ to analyze and
  compare ROC curves''. \emph{BMC Bioinformatics}, \bold{7}, 77.
  DOI: \doi{10.1186/1471-2105-12-77}.
}

\author{Reza Mohammadi \email{a.mohammadi@uva.nl}}

\seealso{ \code{\link[pROC:auc]{pROC::auc()}}, \code{\link[pROC:plot.roc]{pROC::plot.roc()}}, \code{\link[pROC:print.roc]{pROC::print.roc()}}, \code{\link{bdgraph}} }

\examples{
\dontrun{
set.seed( 100 )

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 50, p = 6, size = 7, vis = TRUE )

# Running sampling algorithm based on GGMs 
sample.ggm <- bdgraph( data = data.sim, method = "ggm", iter = 10000 )

# ROC curve for GGM method
roc.ggm <- roc( pred = sample.ggm, actual = data.sim, plot = TRUE )

# Running sampling algorithm based on GCGMs
sample.gcgm <- bdgraph( data = data.sim, method = "gcgm", iter = 10000 )

# ROC curve for GGM and GCGM methods
roc.gcgm <- roc( pred = sample.gcgm, actual = data.sim, plot = TRUE )

ggroc( list( roc.ggm = roc.ggm, roc.gcgm = roc.gcgm ) )
}
}

\keyword{structure learning}

