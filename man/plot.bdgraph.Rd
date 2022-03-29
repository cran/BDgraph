\name{plot.bdgraph}
\alias{plot.bdgraph}

\title{ Plot function for \code{S3} class "\code{bdgraph}" }

\description{ Visualizes structure of the selected graphs which could be a graph with links for which their estimated posterior probabilities are greater than 0.5
or graph with the highest posterior probability. }

\usage{ \method{plot}{bdgraph}( x, cut = 0.5, number.g = NULL, main = NULL, 
           layout = igraph::layout_with_fr, vertex.size = 2, vertex.color = "orange", 
           vertex.frame.color = "orange", vertex.label = NULL, vertex.label.dist = 0.5, 
           vertex.label.color = "blue", edge.color = "lightblue", ... ) }

\arguments{
  \item{x}{ object of \code{S3} class "\code{bdgraph}", from function \code{\link{bdgraph}}. }
  
  \item{cut}{ threshold for including the links in the selected graph based on the estimated posterior probabilities of the links; See the examples. } 

  \item{number.g}{ number of graphs with the highest probabilities. This option works for the case running function \code{bdgraph()} with option \code{save = TRUE}; See the examples. }
  
  \item{main}{ Graphical parameter (see plot).}

  \item{layout}{ vertex placement which is according to \code{R} package \code{\link[igraph]{igraph}}; For different layouts, see \code{\link[igraph]{layout}} of \code{R} package \code{igraph}. }
  \item{vertex.size}{ vertex size which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  \item{vertex.color}{ vertex color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  \item{vertex.frame.color}{ vertex frame color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{vertex.label}{ vertex label. The default vertex labels are the vertex ids. }
  \item{vertex.label.dist}{ vertex label distance which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  \item{vertex.label.color}{ vertex label color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{edge.color}{ edge color which is according to \code{R} package \code{\link[igraph]{igraph}}. }
  
  \item{\dots}{ additional plotting parameters. For the complete list, see \code{\link[igraph]{igraph.plotting}} of \code{R} package \code{igraph}. }
}

\references{
Mohammadi, R. and Wit, E. C. (2019). \pkg{BDgraph}: An \code{R} Package for Bayesian Structure Learning in Graphical Models, \emph{Journal of Statistical Software}, 89(3):1-30, \doi{10.18637/jss.v089.i03} 

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}

Mohammadi, R., Massam, H. and Letac, G. (2021). Accelerating Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Journal of the American Statistical Association}, \doi{10.1080/01621459.2021.1996377} 
  
Dobra, A. and Mohammadi, R. (2018). Loglinear Model Selection and Human Mobility, \emph{Annals of Applied Statistics}, 12(2):815-845, \doi{10.1214/18-AOAS1164}

Mohammadi, A. et al (2017). Bayesian modelling of Dupuytren disease by using Gaussian copula graphical models, \emph{Journal of the Royal Statistical Society: Series C}, 66(3):629-645, \doi{10.1111/rssc.12171} 
}

\author{ Reza Mohammadi \email{a.mohammadi@uva.nl} and Ernst Wit }

\seealso{ \code{\link{bdgraph}}, \code{\link{bdgraph.mpl}} }

\examples{
\dontrun{
set.seed( 100 )

# Generating multivariate normal data from a 'random' graph
data.sim <- bdgraph.sim( n = 100, p = 15, graph = "random", prob = 0.2, vis = TRUE )
   
bdgraph.obj <- bdgraph( data = data.sim )
   
plot( bdgraph.obj )
   
bdgraph.obj <- bdgraph( data = data.sim, save = TRUE )
   
plot( bdgraph.obj, cut = 0.5 )

plot( bdgraph.obj, number.g = 4 )
}
}

\keyword{hplot}

