\name{select}
\alias{select}

\title{Selecting the best graphs}

\description{
According to output of the BD-MCMC algorithm, this function gives us the graphs with the highest posterior probabilities.  
}

\usage{select(output, vis = FALSE)}

\arguments{
  \item{output}{a list which is the result of the BD-MCMC algorithm from the 'bdgraph' function. }
  
  \item{vis}{logical: if TRUE you will see the plot of best graph. The default is FALSE. }
}
\value{
\item{G}{ Adjacency matrix for the best graph (graph with the highest posterior probability).}
}
\references{
Mohammadi, A. and E. C. Wit (2012). Gaussian graphical model determination based on birth-death 
MCMC inference, arXiv:1210.5371v4. \url{http://arxiv.org/abs/1210.5371v4}
}
\author{Abdolreza Mohammadi and Ernst Wit}

\seealso{\code{\link{bdgraph}}}

\examples{
\dontrun{
  # generating the data from a random graph
  data.sim <- bdgraph.sim(n = 100, p = 8, size = 10, vis = TRUE)
  
  output   <- bdgraph(data.sim, mean = 0, iter = 2000)
  
  select(output)
  }
}
