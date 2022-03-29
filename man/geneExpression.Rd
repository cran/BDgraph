\name{geneExpression}
\alias{geneExpression}

\docType{data}

\title{ Human gene expression dataset }

\description{ 
The dataset contains human gene expression of \eqn{100} transcripts (with unique Illumina TargetID) measured on \eqn{60} unrelated individuals. 
}

\usage{ data( geneExpression ) }

\format{ The format is a matrix with 60 rows (number of individuals) and 100 column (number of transcripts). }

\references{
Bhadra, A. and Mallick, B. K. (2013). Joint High Dimensional Bayesian Variable and Covariance Selection with an Application to eQTL Analysis, \emph{Biometrics}, 69(2):447-457, \doi{10.1111/biom.12021}

Mohammadi, A. and Wit, E. C. (2015). Bayesian Structure Learning in Sparse Gaussian Graphical Models, \emph{Bayesian Analysis}, 10(1):109-138, \doi{10.1214/14-BA889}
}

\examples{
data( geneExpression )

dim( geneExpression )
head( geneExpression )
}

\keyword{ datasets }
