\name{geneExpression}
\alias{geneExpression}
\docType{data}
\title{ Human gene expression dataset }

\description{ 
The dataset contains human gene expression of \eqn{100} transcripts 
(with unique Illumina TargetID) measured on \eqn{60} unrelated individuals. 
}

\usage{ data( geneExpression ) }

\format{ The format is a matrix with 60 rows (number of individuals) and 100 column (number of transcripts). }

\source{ The genotypes of those 60 unrelated individuals are available from the Sanger Institute 
          website at ftp://ftp.sanger.ac.uk/pub/genevar
       }

\examples{
\dontrun{
  data( geneExpression )
  dim( geneExpression )
  head( geneExpression )
  boxplot( geneExpression )
}
}

\keyword{ datasets }
