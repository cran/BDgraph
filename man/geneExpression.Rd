\name{geneExpression}
\alias{geneExpression}
\docType{data}
\title{ Human gene expression dataset }

\description{ 
This data set contains human gene expression of \eqn{p = 100} transcrips 
(with unique Illumina TargetID) measured on \eqn{n = 60} unrelated individuals. 
}

\usage{ data(geneExpression) }

\format{ The format is a matrix with 60 rows (number of individuals) and 100 column (number of transcrips). }

\source{ The genotypes of those 60 unrelated individuals are available from the Sanger Institute 
          website at ftp://ftp.sanger.ac.uk/pub/genevar
       }

%\author{Abdolreza Mohammadi and Ernst Wit}

\examples{
  data( geneExpression )
  dim( geneExpression )
  head( geneExpression )
  boxplot( geneExpression )
}

\keyword{ datasets }