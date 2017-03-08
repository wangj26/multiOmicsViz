multiOmicsViz <- function(sourceOmics,sourceOmicsName,chrome_sourceOmics,
targetOmicsList,targetOmicsName,chrome_targetOmics,fdrThr,outputfile,
nThreads=NULL,legend=TRUE){
    
    outputfile <- paste(outputfile,".png",sep="")
    
    if(class(sourceOmics)=="SummarizedExperiment"){
      sourceOmics <- assays(sourceOmics,n=1)
    }else{
      if(!is.matrix(sourceOmics) && !is.data.frame(sourceOmics)){
        stop("Source omics data (e.g. CNA data) should be a R matrix, 
        data.frame or SummarizedExperiment object.")
      }
    }
    
    if(!is.character(sourceOmicsName)){
      stop("sourceOmicsName is the name of the source omics data, which 
      should be a R character object.") 
    }
    
    if(!is.list(targetOmicsList)){
      stop("Multipl target omics data (e.g. mRNA or protein data) 
      should be saved in the list object.")
    }
    
    if(length(targetOmicsList)>5){
      stop("targetOmicsList can only contain at most five omics data.")
    }
    
    for(i in seq_len(length(targetOmicsList))){
      if(class(targetOmicsList[[i]])=="SummarizedExperiment"){
        targetOmicsList[[i]] <- assays(targetOmicsList[[i]],n=1)
      }else{
        if(!is.matrix(targetOmicsList[[i]]) && 
        !is.data.frame(targetOmicsList[[i]])){
          stop("Each of all target omics data in the list (e.g. mRNA or 
          protein data) should be a R matrix, data.frame or 
          SummarizedExperiment object.")
        }
      }
    }
    
    if(!is.character(targetOmicsName)){
      stop("targetOmicsName should be a R vector object 
      containing the name for each omics data in the targetOmicsList.")
    }
    
    if(length(targetOmicsList)!=length(targetOmicsName)){
      stop("targetOmicsName should have the same length 
      with targetOmicsList.")
    }   
    
    if(length(targetOmicsList)>1 && is.null(nThreads)){
      stop("Please input nThreads for the parallel computing.")
    }
    
    if(length(targetOmicsList)>1 && !is.null(nThreads)){
    
      if(nThreads>length(targetOmicsList)){
        stop("nThreads should be at most the length of targetOmicsList.")
      }
    }
    
        
    ########find the intersect genes among all omics data######
    intG <- c()
    for(i in seq_len(length(targetOmicsList))){
      if(i==1){
        intG <- rownames(targetOmicsList[[i]])
      }else{
        intG <- intersect(intG,rownames(targetOmicsList[[i]]))
      }
    }
    
    if(length(intG)==0){
      stop("The ID types of all omics data in the targetOmicsList 
      should be the same.")
    }
    
    #####process chrome location###
    
    chromeList <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
    "14","15","16","17","18","19","20","21","22","X","Y","All")
    
    x <- setdiff(chrome_sourceOmics,chromeList)
    if(length(x)>0){
      stop('The input chrome infomation for source omics data contains the 
      invalid information. Please only input the chromosome names from 
      the following list: "1","2","3","4","5","6","7","8","9","10","11","12",
      "13","14","15","16","17","18","19","20","21","22","X","Y" and "All".')
    }
    
    x <- setdiff(chrome_targetOmics,chromeList)
    if(length(x)>0){
      stop('The input chrome infomation for target omics data contains the 
      invalid information. Please only input the chromosome names from 
      the following list: "1","2","3","4","5","6","7","8","9","10","11","12",
      "13","14","15","16","17","18","19","20","21","22","X","Y" and "All".')
    }

    if((length(chrome_sourceOmics)>1 && which(chrome_sourceOmics=="All")>1) 
    || (length(chrome_sourceOmics)==1 && chrome_sourceOmics=="All")){
      chrome_sourceOmics <- "All"
    }
    
    if((length(chrome_targetOmics)>1 && which(chrome_targetOmics=="All")>1) 
    || (length(chrome_targetOmics)==1 && chrome_targetOmics=="All")){
      chrome_targetOmics <- "All"
    }
    
    if(chrome_sourceOmics=="All"){
      chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
      "12","13","14","15","16","17","18","19","20","21","22","X","Y")
    }
      
    if(chrome_targetOmics=="All"){
      chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
      "12","13","14","15","16","17","18","19","20","21","22","X","Y")
    }
      
    #######Extract sub list#########
    genelocate <- datacache$genelocate
    
    genelocate_sourceOmics <- genelocate[genelocate[,2] %in%
chrome_sourceOmics,]
    genelocate_targetOmics <- genelocate[genelocate[,2] %in%
chrome_targetOmics,]
  
    intG <- intersect(intG,genelocate_targetOmics[,1])
    
    if(length(intG)==0){
      stop("The ID types for all omics data in the targetOmicsList should be 
      gene symbol or all genes in the target omics data are not in the 
      selected chromosomal location chrome_targetOmics.")
    }
    
    for(i in seq_len(length(targetOmicsList))){
      targetOmicsList[[i]] <- targetOmicsList[[i]][intG,]
    }
        
    source_gene <- rownames(sourceOmics)
    source_gene_locate <-
intersect(unique(genelocate_sourceOmics[,1]),source_gene)
    if(length(source_gene_locate)==0){
      stop("The ID type in the source omics data should be gene symbol or all 
      genes in the source omics data are not in the selected chromosomal 
      location chrome_sourceOmics.")
    }
    source_gene <- sourceOmics[source_gene_locate,]
    
    genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[,1] 
    %in% source_gene_locate,]
    genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[,1] 
    %in% intG,]
  
    
    ###Calculate the correlation between cna and other omics data######
    cat("Identify the significant correlations...\n")
    if(length(targetOmicsList)==1){
      resultList <-
calculateCorForTwoMatrices(source_gene,targetOmicsList[[1]],fdrThr)
    }else{
      cl <- makeCluster(nThreads)
      registerDoParallel(cl)
      resultList <- list()
      resultList <- foreach(i=seq_len(length(targetOmicsList)), 
      .packages="multiOmicsViz") %dopar% {
        corrArray <-
calculateCorForTwoMatrices(source_gene,targetOmicsList[[i]],fdrThr)
        return(corrArray)
      }
      stopCluster(cl)
    }
        
    
    ##Calculate the location of genes in the heatmap
    chromLength <- datacache$chromLength
    
    re <-
.calculateChromLength(chromLength,chrome_sourceOmics,genelocate_sourceOmics)
    genelocate_sourceOmics <- re$genelocate
    chromLength_sourceOmics <- re$chromLength
    
    re <-
.calculateChromLength(chromLength,chrome_targetOmics,genelocate_targetOmics)
    genelocate_targetOmics <- re$genelocate
    chromLength_targetOmics <- re$chromLength
    
    ##########Plot Figure############
    
    cat("Plot figure...\n")
    
    if(length(targetOmicsList)==1){
      png(outputfile,height=480*5,width=480*5,res=300)
      .plotHeatMap(resultList,genelocate_sourceOmics,chromLength_sourceOmics,
      genelocate_targetOmics,chromLength_targetOmics,sourceOmicsName,
      targetOmicsName,dim=1)
      if(legend==TRUE){
        legend("topleft",c("positive correlation","negative correlation"),
        col=c("red","green"),pch=19)
      }
    }else{
      png(outputfile,height=480*8,width=480*5*length(targetOmicsList),res=300)
      layout(matrix(c(1:(2*length(targetOmicsList))),length(targetOmicsList),
      2,byrow=TRUE),heights=c(2,1))
      for(i in seq_len(length(resultList))){
       
.plotHeatMap(resultList[[i]],genelocate_sourceOmics,chromLength_sourceOmics,
        genelocate_targetOmics,chromLength_targetOmics,sourceOmicsName,
        targetOmicsName[i],dim=2)
      }
      if(legend==TRUE){
        legend("topleft",c("positive correlation","negative correlation"),
        col=c("red","green"),pch=19)
      }
     
.plotSummaryBar(resultList,chromLength_sourceOmics,genelocate_sourceOmics,
      sourceOmicsName)
      if(legend==TRUE){
        legend("topleft",c("specific correlation","common correlation"),
        col=c("blue","black"),lty=1,lwd=3)
      }
    }
    dev.off()
}

.calculateChromLength <- function(chromLength,selectedChrom,genelocate){
    chromLength <- chromLength[chromLength[,1] %in% selectedChrom,,drop=FALSE]
    
    if(length(selectedChrom)==1){
      x <- 0
    }else{
      x <- c(0,chromLength[1:(nrow(chromLength)-1),2])
    }
    chromLength[,3] <- cumsum(as.numeric(x))
    chromLength[,4] <- cumsum(as.numeric(chromLength[,2]))
    
    genelocate <- cbind(genelocate,0,0)
    
    colnames(genelocate)[5:6] <- c("finalstart","finalend")
        
    for(i in c(1:nrow(genelocate))){
        chr <- genelocate[i,2]
        s <- genelocate[i,3]
        e <- genelocate[i,4]
        cs <- chromLength[chromLength[,1]==chr,3]
        genelocate[i,5] <- s+cs
        genelocate[i,6] <- e+cs
    }
    re <- list(chromLength=chromLength,genelocate=genelocate)
    return(re)
}

.plotHeatMap <- function(corrArray,genelocate_sourceOmics,
chromLength_sourceOmics,genelocate_targetOmics,chromLength_targetOmics,
sourceOmicsName,targetOmicsName,dim=1){

    allChromlen_sourceOmics <-
    chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]
    allChromlen_targetOmics <-
    chromLength_targetOmics[nrow(chromLength_targetOmics),4]
    
    if(dim==1){
      par(mar=c(4,4,4,0))
    }else{
      par(mar=c(0,4,4,0))
    }

    p <- which(corrArray!=0,arr.ind=TRUE)
    allcnagene <- rownames(corrArray)
    allovgene <- colnames(corrArray)

    la <- 1
    for(i in c(1:nrow(p))){
       
        cnag <- allcnagene[p[i,1]]
        ovg <- allovgene[p[i,2]]
        cnagp <- genelocate_sourceOmics[genelocate_sourceOmics[,1]==cnag,5]
        ovgp <- genelocate_targetOmics[genelocate_targetOmics[,1]==ovg,5]
    
        if(length(cnagp)==0 || length(ovgp)==0){
            next
        }
    
        cov <- corrArray[cnag,ovg]
        color <- ifelse(cov>0,"red","green")
       
        if(la==1){
            if(dim==1){
              plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", targetOmicsName,
              " correlation",sep=""), xlim=c(0,allChromlen_sourceOmics),
              ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
              frame.plot=FALSE,xlab=paste(sourceOmicsName,
              " chromosomal location",sep=""),ylab=paste(targetOmicsName,
              " chromosomal location",sep=""),pch=20,col=color,cex=0.2)
              axis(side=1,at=(chromLength_sourceOmics[,4]-
chromLength_sourceOmics[,2]/2),
              labels=chromLength_sourceOmics[,1])
            }else{
              plot(cnagp,ovgp,main=paste(sourceOmicsName,"-", 
              targetOmicsName," correlation",sep=""),
xlim=c(0,allChromlen_sourceOmics),
              ylim=c(0,allChromlen_targetOmics),xaxt="n",yaxt="n",
              frame.plot=FALSE,ylab=paste(targetOmicsName," chromosomal
location",sep=""),
              xlab="",pch=20,col=color,cex=0.2)
            }
            axis(side=2,at=(chromLength_targetOmics[,4]-
chromLength_targetOmics[,2]/2),
            labels=chromLength_targetOmics[,1])
           
abline(h=c(0,chromLength_targetOmics[,4]),v=c(0,chromLength_sourceOmics[,4]),
            col="gray",lty=3)
            la <- la+1
        }else{
            for(u in seq_len(length(cnagp))){
                for(v in seq_len(length(ovgp))){
                    points(cnagp[u],ovgp[v],pch=20,col=color,cex=0.2)
                }
            }
        }
    }
}

.plotSummaryBar <- function(resultList,chromLength_sourceOmics,
genelocate_sourceOmics,sourceOmicsName){
  #######Identify the common pairs for all omics data
    allChromlen <- chromLength_sourceOmics[nrow(chromLength_sourceOmics),4]

    sumA <- abs(resultList[[1]])
    for(i in seq(from=2,to=length(resultList),by=1)){
      sumA <- sumA+abs(resultList[[i]])
    }
    ov <- sumA
    ov[ov!=length(resultList)] <- 0
    ov <- apply(ov,1,sum)
    ov <- (-ov)/length(resultList)
    
    maxO <- min(ov)
    
    ######Identify specific pairs for all omics data
    spe <- list()
    sumA[sumA>1] <- 0
    maxP <- 0
    for(i in seq_len(length(resultList))){
      y <- abs(resultList[[i]])*sumA
      spe[[i]] <- apply(y,1,sum)
      maxP <- max(maxP,max(spe[[i]]))
    }
    
    for(i in seq_len(length(spe))){
    
      par(mar=c(4,4,0,0))

      plot(0,0,xlim=c(0,allChromlen),ylim=c(maxO,maxP),type="n",xaxt="n",
      frame.plot=FALSE,xlab=paste(sourceOmicsName," chromosomal
location",sep=""),
      ylab="Number of significant correlations")
      
      axis(side=1,at=(chromLength_sourceOmics[,4]-
chromLength_sourceOmics[,2]/2),
      labels=chromLength_sourceOmics[,1])
      axis(side=2,at=NULL,labels=TRUE)
      abline(v=c(0,chromLength_sourceOmics[,4]),col="gray",lty=3)
      
      spe_o <- spe[[i]]
      
      for(i in seq_len(length(spe_o))){
        gg <- names(spe_o)[i]
        gg_p <- genelocate_sourceOmics[genelocate_sourceOmics[,1]==gg,5]
        if(length(gg_p)==0){
            next
        }
        up <- spe_o[gg]
        do <- ov[gg]
        for(j in seq_len(length(gg_p))){
            points(gg_p[j],up,cex=0.2,type="h",col="blue")
            points(gg_p[j],do,cex=0.2,type="h",col="black")
        }
      }
    }
}