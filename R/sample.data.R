sample.data <- function(Data, cleanData, numbeRow=300){
    
    sampled <- list()
    ##set.seed(seed)
    ##^Set the seed outside,
    ##https://support.bioconductor.org/p/110439/
    
    ## Sampling raw Data
    ## Locus2gene
    locusInfo <- as.data.frame(SummarizedExperiment::rowRanges(Data$rawDnam))
    lociSampled1 <- sample(rownames(cleanData$locus2gene), size=numbeRow*2, 
                           replace=FALSE)
    randomLoci <- locusInfo[!locusInfo$probeID %in% lociSampled1, "probeID"]
    lociSampled2 <- sample(randomLoci, size=numbeRow*2, replace=FALSE)
    sampleDnam <- Data$rawDnam[c(lociSampled1, lociSampled2), ]
    locusInfoSampled <- as.data.frame(SummarizedExperiment::rowRanges(sampleDnam))  
    geneSampled <- unique(locusInfoSampled$Gene_Symbol)  
    genExprRaw <- Data$genExpr[rownames(Data$genExpr) %in% geneSampled, ]
    rawData <- list(genExpr=genExprRaw, genExprSampleInfo=Data$genExprSampleInfo, 
                    rawDnam=sampleDnam, clinical=Data$clinical)
    sampled[["rawData"]] <- rawData 

    ## sampling cleaned data
    l2g <- cleanData$locus2gene[intersect(rownames(cleanData$locus2gene), 
                                          c(lociSampled1, lociSampled2)), ]
    dnam <- cleanData$dnam[rownames(cleanData$dnam) %in% rownames(l2g), ]
    genExpr <- cleanData$genExpr[rownames(cleanData$genExpr) %in% 
                                       unique(l2g$Gene_Symbol), ]
    
    cleaned <- list(dnam=dnam, genExpr=genExpr, survival=cleanData$survival, 
                    locus2gene=l2g)

    sampled[["cleaned"]] <- cleaned
    return(sampled)
}
