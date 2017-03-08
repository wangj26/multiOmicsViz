\name{multiOmicsViz}
\alias{multiOmicsViz}

\title{
	Plot the effect of one omics data on other omics data along the chromosome
}

\description{
	Calculate the spearman correlation between the source omics data and other 
	target omics data, identify the significant correlations and plot the significant 
	correlations on the heat map in which the x-axis and y-axis are ordered by the 
	chromosomal location.
}

\usage{
multiOmicsViz(sourceOmics,sourceOmicsName,chrome_sourceOmics,targetOmicsList,
targetOmicsName,chrome_targetOmics,fdrThr,outputfile,nThreads=NULL,legend=TRUE)
}

\arguments{
  \item{sourceOmics}{
	A R matrix, data.frame or SummarizedExperiment object containing the omics data.
	The data should contain the row names representing the genes and column names 
	representing the samples.
  }
  \item{sourceOmicsName}{
	The name of the source omics data (e.g. CNA).
  }
  \item{chrome_sourceOmics}{
	The \code{multiOmicsViz} function will extract the genes in the selected 
	chromosome(s) from genes in the source omics data and then identify and visualize the 
	significant correlations based on the selected genes. \code{chrome_sourceOmics}
	can be one character containing the chromosome name (e.g. "1"), a R vector object
	containing multiple chromosomes (e.g. c("1","2","3")) or "All" representing all 
	24 chromosomes.
  }
  \item{targetOmicsList}{
	A R list object containing at most 5 target omics data. 
	Each omics data in the list should be a R matrix, data.frame or SummarizedExperiment 
	object and contain the row names representing the genes and column names 
	representing the samples. There should have multiple overlapping genes among all target 
	omics data and at least 6 overlapping sample between source omics data and each target 
	omics data.
  }
  \item{targetOmicsName}{
  A R vector object containing the name of all target omics data stored in 
  the \code{targetOmicsList}.
  }
  \item{chrome_targetOmics}{
	The \code{multiOmicsViz} function will extract the genes in the selected 
	chromosome(s) from the overlapping genes among all target omics data and 
	then identify and visualize the significant correlations based on the selected genes. 
	\code{chrome_sourceOmics} can be one character containing the chromosome name 
	(e.g. "1"), a R vector object containing multiple chromosomes (e.g. c("1","2","3")) 
	or "All" representing all 24 chromosomes.
  }
  \item{fdrThr}{
  The FDR threshold for identifying the significant correlations.
  }
  \item{outputfile}{
	The output file name. 
  }
  \item{nThreads}{
	If \code{targetOmicsList} contains 2 or 3 omics data, multiOmicsViz will use
	the parallel computing method to calculate the significant correlations between
	the source omics data and each of target omics data. \code{nThreads} is the number
	of cores used for the parallel computing.
  }
  \item{legend}{
	If \code{legend} is TRUE, the output heat map will contain the legend.
  }
}

\value{
	If the \code{targetOmicsList} contains one target omics data, the \code{multiOmicsViz}
	function will plot a heat map in which x-axis represents the genes in the source omics 
	data, y-axis represents the genes in the target omics data, x-axis and y-axis are 
	ordered by chromosomal location, each point represents a significant correlation, red 
	color represents the significant positive correlation and blue color represents the 
	significant negative correlation.
	If the \code{targetOmicsList} contains multiple target omics data, the 
	\code{multiOmicsViz} function will not only plot multiple heat maps for each target 
	omics data but also plot mutliple bar charts in which blue bars represent the number 
	of specific significant correlations for the target omics data and black bars 
	represents the number of common significant correlations among all target omics data.
}

\author{
	Jing Wang
}

\examples{
	sourceOmics <- system.file("extdata","sourceOmics.txt",package="multiOmicsViz")
	sourceOmics <- read.table(sourceOmics,header=TRUE,sep="\t",stringsAsFactors=FALSE,
	check.names=FALSE)
	
	targetOmics1 <- system.file("extdata","targetOmics.txt",package="multiOmicsViz")
	targetOmics1 <- read.table(targetOmics1,header=TRUE,sep="\t",stringsAsFactors=FALSE,
	check.names=FALSE)
	
	targetOmicsList <- list()
	targetOmicsList[[1]] <- targetOmics1
	
	outputfile <- paste(tempdir(),"/heatmap",sep="")
	multiOmicsViz(sourceOmics,"CNA","20",targetOmicsList,
"mRNA","All",0.001,outputfile)
}

\keyword{ methods }
