\name{transfer}
\alias{transfer}

\title{ transfer for discrete data }

\description{ Transfers discrete data, by counting the duplicated rows.  }

\usage{ transfer( r_data ) }

\arguments{
  \item{r_data}{ An (\eqn{n \times p}{n x p}) matrix or a \code{data.frame} corresponding to the data (\eqn{n} is the sample size and \eqn{p} is the number of variables). }
}

\value{
\item{data}{ An (\eqn{n \times p+1}) matrix of transferred data, in which the last column is the frequency of duplicated rows. }
}

%\references{ }

\author{ Abdolreza Mohammadi }

\seealso{ \code{\link{bdgraph.sim}}, \code{\link{bdgraph.mpl}} and \code{\link{bdgraph}} }

\examples{
\dontrun{
# Generating multivariate binary data from a 'random' graph
data.sim <- bdgraph.sim( n = 12, p = 4, size = 4, type = "discrete", cut = 2 )
r_data     <- data.sim $ data
r_data
   
# Transfer the data
transfer( r_data )
  
}
}
