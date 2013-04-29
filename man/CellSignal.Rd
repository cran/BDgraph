\name{CellSignal}
\alias{CellSignal}
\docType{data}
\title{A flow cytometry dataset}
\description{
This data set contains flow cytometry of \eqn{p = 11}{p = 11} proteins measured on \eqn{n = 11672}{n = 11672} cells.
}
\usage{data(CellSignal)}
\format{
  The format is a list with 11672 rows (number of cells) and 11 column (number of proteins).
}
\source{
It is publicly available at http://www.sciencemag.org/content/308/5721/523/suppl/DC1
}
\author{
Abdolreza Mohammadi and Ernst Wit
}
\seealso{
\code{\link{bdgraph}}
}
\examples{
\dontrun{
  data(CellSignal)
  head(CellSignal)
  }
}

\keyword{datasets}