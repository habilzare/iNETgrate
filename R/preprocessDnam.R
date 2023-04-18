preprocessDnam <- function(rawDnam, rawDnamSampleInfo=NULL, savePath, 
                            annLib="Auto", verbose=0){
    ## This function preprocess the dnam data to remove probes with >50% of
    ## missing Beta values (NAs), SNP-enriched probes, non-CpG probes and
    ## imputes probes with few missing Beta values.
    ## Input:
    ## rawDnam: Summarized experiment object obtained similar to TCGABiolinks prepare
    ## function. Here assay(rawDnam) is a matrix of Beta values.
    ## rawDnamSampleInfo: A dataframe of rawDnam sample information
    ## annLib: Annotation library to be used for identifying SNPs and non-Cpg
    ## probes. Default is set to "IlluminaHumanMethylation450kanno.ilmn12.hg19".
    ## Output:
    ## A list of containing:Beta matrix, locus2gene dataframe, sample
    ## information dataframe, sample2patient vector 

    result <- list()
    starTime <- Sys.time()

    message.if(me=paste("Preprocessing dnam started at:", 
                      format(starTime, usetz=TRUE),"\n"), verbose=verbose-3)
    message.if(paste("preprocessDnam results will be saved at:", 
                   savePath, "\n"), verbose=verbose-1)

    if(annLib=="Auto")
        annLib <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"

    ann450k <- minfi::getAnnotation(annLib)
    ## QC:
    ## c1 <- class(rawDnam)
    ## If rawDnam is matrix, we prepare a SummarizedExperiment object
    if(inherits(rawDnam, "matrix")){
        message.if(me=paste("rawDnam is the beta matrix!!\n"), verbose=verbose)
        subAnn <- ann450k[ ,c("chr", "pos", "Name", "UCSC_RefGene_Name")]
        ##^Later this may need a fix. Not sure if all the annotation libraries
        ##have the above mentioned column names.
        colnames(subAnn) <- c("Chromosome", "Genomic_Coordinate", "probeID","Gene_Symbol")
        common <- intersect(rownames(rawDnam), rownames(subAnn))
        subAnn <- subAnn[common, ]
        gr1 <- GenomicRanges::makeGRangesFromDataFrame(subAnn, keep.extra.columns=TRUE,
               start.field="Genomic_Coordinate", end.field="Genomic_Coordinate")
        se <- SummarizedExperiment(rawDnam[common,], rowRanges=gr1)
        if(!is.null(rawDnamSampleInfo)){
            se <- SummarizedExperiment(rawDnam[common,], rowRanges=gr1,
                                       colData=rawDnamSampleInfo[colnames(rawDnam), ])
        } 
        rawDnam <- se
    }

    ## Remove rs or sex probe IDs  (probes that donot have assigned chromosome info)
    message.if(me="Step 1: Removing loci that donot have assigned chromosome...",
               verbose=verbose-1)
    chrNA <- as.character(GenomicRanges::seqnames(rawDnam)) %in% c("chrNA", "chrX", "chrY") 
    subData1 <- subset(rawDnam, subset=!chrNA)
    message.if(me=paste(sum(chrNA),"loci with no chromosome assignment were removed.\n"), 
               verbose=verbose-1)


    ## Remove probes that have more than 50% NAs
    message.if(me="Step 2: Removing loci that are NA in more than half of the samples...",
               verbose=verbose-1)
    naNum <- rowSums(is.na(SummarizedExperiment::assay(subData1)))
    naPerSample <- colSums(is.na(SummarizedExperiment::assay(subData1)))
    plotPath <- file.path(savePath, "plots")
    dir.create(plotPath, showWarnings=FALSE)
    png(file.path(plotPath,"nas_per_loci.png"), height=480, width=2*480)
    par(mfrow=c(1,2), cex.axis=1.2, cex.lab=1.5)
    plot(x=sort(naNum), y=seq_along(naNum), xlab="# of NA per locus",
         ylab="# of loci")
    plot(x=sort(naPerSample), y=seq_along(naPerSample), 
         xlab="# of NA per sample", ylab="# of samples")
    dev.off()
    tooManyNAs <- naNum > ncol(subData1)/2
    subData2 <- subset(subData1, subset=!tooManyNAs)
    message.if(me=paste(sum(tooManyNAs),"loci with more than 50% NAs were removed.\n"), 
               verbose=verbose-1)

    ## Remove non-CpG (CH) probes
    message.if(me="Step 3: Removing non-CpG Loci...", verbose=verbose-1)
    chProbes <- rownames(ann450k)[substring(rownames(ann450k), 1,2)!="cg"]
    toKeep1 <- !rownames(SummarizedExperiment::rowData(subData2)) %in% chProbes
    names(toKeep1) <- rownames(SummarizedExperiment::rowData(subData2))
    subData3 <- subset(subData2, subset=toKeep1)
    message.if(me=paste(sum(!toKeep1),"non-CpG loci were removed.\n"),
               verbose=verbose-1)

    ## Remove SNP enriched probes
    message.if(me="Step 4: Removing SNP enriched Loci...", verbose=verbose-1)
    snpCols <- c("Probe_rs","Probe_maf","CpG_rs","CpG_maf", "SBE_rs", "SBE_maf")
    sbeSNP <- rownames(ann450k)[which(ann450k[,"SBE_maf"] != "NA")]
    cpgSNP <- rownames(ann450k)[which(ann450k[,"CpG_maf"] != "NA")]
    removeSNPs <- union(sbeSNP, cpgSNP)
    toKeep2 <- !rownames(SummarizedExperiment::rowData(subData3)) %in% removeSNPs
    names(toKeep2) <- rownames(SummarizedExperiment::rowData(subData3))
    subData4 <- subset(subData3, subset=toKeep2)
    message.if(me=paste(sum(!toKeep2),"loci with SNP enrichment were removed.\n"),
               verbose-1)

    ## Impute probes with <50% NAs.
    message.if(me="Step 5: Imputing NAs in beta matrix with mean beta values ...", 
               verbose=verbose-1)
    if(sum(is.na(SummarizedExperiment::assay(subData4)))/length(SummarizedExperiment::assay(subData4))>0.01){
        stop("More than 1% NA in betaMatrix!")
    } else { ## Replacing NAs with the average of beta in that loci
        nas <- which(is.na(SummarizedExperiment::assay(subData4)))
        meanMatrix <- matrix(rowMeans(SummarizedExperiment::assay(subData4), na.rm=TRUE),
                             ncol=ncol(SummarizedExperiment::assay(subData4)),
                             nrow=nrow(SummarizedExperiment::assay(subData4)), byrow=FALSE)
        message.if(me=paste(length(nas), "NAs in",
                    sum(is.na(rowMeans(SummarizedExperiment::assay(subData4), na.rm=FALSE))), 
                    "loci, which will be replaced by the mean beta value of the
                    corresponding loci.\n"), verbose=verbose-1)
        SummarizedExperiment::assay(subData4)[nas] <- meanMatrix[nas]
    }
    ## Double check NAs:
    if(any(is.na(assay(subData4))))
        stop("betaMatrix is not expected to have NAs!")

    dnam <- SummarizedExperiment::assay(subData4)
    result[["dnam"]] <- dnam
    sampleInfo <- SummarizedExperiment::colData(subData4)
    result[["sampleInfo"]] <- sampleInfo
    ## s2pDnam <- structure(sampleInfo[,"patient"], names=sampleInfo[,"barcode"])
    ## result[["s2p"]] <- s2pDnam
    
    ## Prepare locus2gene dataframe by adding RefGene_Group
    message.if(me="Preparing locus2gene...", verbose=verbose-1)
    locusInfo <- as.data.frame(SummarizedExperiment::rowRanges(subData4))
    ## locus2gene <- locusInfo[FALSE,c("probeID","Gene_Symbol","seqnames", "start")]
    if(nrow(locusInfo[locusInfo$width != 1,]) !=0){
        warning("Genomic coordinates not defined. Check locus2gene ouput!")
    } else{
        locus2gene <- locusInfo[, c("probeID","Gene_Symbol", "seqnames", "start")]
    }
    colnames(locus2gene) <- c("probeID","Gene_Symbol","Chromosome",
                              "Genomic_Coordinate")
    rownames(locus2gene) <- locus2gene$probeID
    result[["locus2gene"]] <- locus2gene

    timeTaken <- Sys.time()-starTime
    message.if(me=paste("Processing DNA methylation data took:", timeTaken, 
                      attr(timeTaken, "units"), "\n"), verbose=verbose-3)

    return(result)
}
