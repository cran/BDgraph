\name{CellSignal}
\alias{CellSignal}
\docType{data}
\title{A flow cytometry dataset}

\description{ This dataset contains flow cytometry of \eqn{11} proteins measured on \eqn{11672} red blood cells.}

\usage{ data( CellSignal ) }

\format{ The format is a matrix with 11672 rows (number of cells) and 11 column (number of proteins).}

\source{Data are publicly available at http://www.sciencemag.org/content/308/5721/523/suppl/DC1}

%\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
  data( CellSignal )
  head( CellSignal )
  boxplot( CellSignal )
}

\keyword{datasets}
