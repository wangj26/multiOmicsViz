calculateCorForTwoMatrices <- function(matrix1,matrix2,fdr){
	####the row of two matrices should be gene ids while col should be samples
	###two matrices should have the same or similar sample names
	
	ov_sample <- intersect(colnames(matrix1),colnames(matrix2))
	
	if(length(ov_sample)<6){
		stop("The number of overlapping samples of two matrices is less than 6 and thus can not perform correlation analysis. Please check the data.")
	}
	
	matrix1 <- matrix1[,ov_sample]
	matrix2 <- matrix2[,ov_sample]
	
	suppressWarnings(corrArray <- cor(t(matrix1),t(matrix2),method="spearman"))
	corrArray[is.na(corrArray)] <- 0

	n <- t(!is.na(t(matrix1))) %*% (!is.na(t(matrix2)))
  suppressWarnings(t <- (corrArray * sqrt(n - 2))/sqrt(1 - corrArray^2))
  suppressWarnings(corrP <- 2 * (1 - pt(abs(t), (n - 2))))
    
  corrP[is.na(corrP)] <- 0
  corrP[corrP>1] <- 1
  
  corrSig <- apply(corrP,1,.fdrFilter,fdr)
  corrSig <- t(corrSig)
  
  corrArray <- sign(corrArray)*corrSig
  return(corrArray)	
}

.fdrFilter <- function(vector,netFDRThr){
    adjp <- p.adjust(vector,method="BH")
    vector[names(adjp)] <- adjp
    sigP <- ifelse(vector<=netFDRThr,1,0)
    return(sigP)
}

