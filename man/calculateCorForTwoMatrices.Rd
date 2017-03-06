\name{calculateCorForTwoMatrices}
\alias{calculateCorForTwoMatrices}

\title{
	Identify the significant correlations between two matrices.
}

\description{
	The calculateCorForTwoMatrices function uses the spearman correlation to 
	identify the significant correlations between two matrices.
}

\usage{
calculateCorForTwoMatrices(matrix1,matrix2,fdr)
}

\arguments{
  \item{matrix1}{
 	A R matrix, data.frame or SummarizedExperiment object containing 
	the numeric values.
  }
  \item{matrix2}{
 	A R matrix, data.frame or SummarizedExperiment object containing 
	the numeric values. \code{matrix2} should have at least 6 overlapping
	samples with \code{matrix1}.
  }
  \item{fdr}{
  	The FDR threshold for identifying significant correlations.
  }
}

\author{
	Jing Wang
}

\value{
	This function will return a R matrix object containing significant 
	correlations. "1" represents the significant positive correlation, "-1"
	represents the significant negative correlation and "0" represents no 
	significant correlation.
}

\examples{
	matrix1 <- system.file("extdata","sourceOmics.txt",package="multiOmicsViz")
	matrix1 <- read.table(matrix1,header=TRUE,sep="\t",stringsAsFactors=FALSE,
	check.names=FALSE)
	
	matrix2 <- system.file("extdata","targetOmics.txt",package="multiOmicsViz")
	matrix2 <- read.table(matrix2,header=TRUE,sep="\t",stringsAsFactors=FALSE,
	check.names=FALSE)
	
	sig <- calculateCorForTwoMatrices(matrix1=matrix1,
	matrix2=matrix2,fdr=0.01)
}

\keyword{ methods }
