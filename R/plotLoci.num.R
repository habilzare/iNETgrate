plotLoci.num <- function(locus2gene, genesColName="Gene_Symbol", 
                         lociColName="probeID", selectedLoci="Auto", plotFile, 
                         doAddTitle=TRUE, xlab1="Number of loci per gene", 
                         ylab1="Cumulative probability", verbose=0){
    ## Plots the number of loci per gene.
    ## selectedLoci: A character vector of probe IDs, e.g., c("cg00000029", "cg18799241").
    ## Habil changeed locus2gene to locusGene on 2019-06-21.
    ## locus2gene: A matrix with at least 2 columns:
    ## probeID (e.g., cg00000029), and "Gene_Symbol" (symbols e.g., TSPY4).
    ##^If set to "hg19", the required information will be automatically inferred.
    result <- list()
    
    ## selectedLoci:
    ## QC:
    if(length(selectedLoci)==0)
        stop("No selectedLoci!")
    createdLocusGene <- create.locusGene(locus2gene, genesColName=genesColName,
                                         lociColName=lociColName,
                                         verbose=verbose)
    locus2oneGene <- createdLocusGene$locus2oneGene
    result[["noGenesAll"]] <- createdLocusGene$noGenes
    if((selectedLoci=="Auto")[1]){
        selectedLoci <- unique(locus2oneGene[,"probeID"])
        result[["noGenes"]] <- result$noGenesAll
    } else {
        result[["noGenes"]] <- intersect(result$noGenesAll, selectedLoci)
    }

    ## Exluding the loci we are not interested in:
    locus2oneGeneSubset <- locus2oneGene[locus2oneGene[,"probeID"] %in% selectedLoci, , drop=FALSE]
    locus2oneGeneSubset <- locus2oneGeneSubset[!locus2oneGeneSubset[,"probeID"] %in%
                                               result$noGenes, , drop=FALSE]
    l2g <- sort(table(locus2oneGeneSubset[, "Gene_Symbol"]),decreasing=FALSE)
    ## Preparation for the plots:
    x1 <- as.numeric(l2g) ## Number of loci per gene
    y1 <- seq_along(l2g)/length(l2g) ## Cumulative probability
    title1 <- paste("95% of genes have less than",x1[0.95*length(x1)]+1,"loci")
    title2 <- paste(length(result$noGenes),"loci are not related to any gene.")
    message.if(cat(title1, "\n"), verbose=verbose-2)
    message.if(cat(title2, "\n"), verbose=verbose-2)

    plotFile2 <- gsub(plotFile, pattern="\\.png$", replacement="_95.png")
    ##png(plotFile, height=480, width=2*480, res=72)
    ##par(mfrow=c(1,2))
    png(plotFile, height=2*480, width=2*480, res=2.5*72)
    plot(x=x1,y=y1, xlab=xlab1, ylab=ylab1, cex.lab=1.4)
    if(doAddTitle & length(result$noGenes)>0) ## The information is available.
        title(title2)
    dev.off()
    ## The zoomed plot:
    aty <- c(0, 0.2, 0.4, 0.6, 0.8, 0.95)
    png(plotFile2, height=2*480, width=2*480, res=2.5*72)
    plot(x=x1[seq_len(0.95*length(x1))],y=y1[seq_len(0.95*length(y1))], xlab=xlab1, ylab=ylab1,
         cex.lab=1.4)
    axis(2, at=aty, labels=aty)
    if(doAddTitle)
        title(title1)
    dev.off()
    message.if(cat("The loci plots was saved at:", 
                   dirname(plotFile), "\n"), verbose=verbose-1)
    result[["l2g"]] <- l2g
    return(result)
}
