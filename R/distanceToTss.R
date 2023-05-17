distanceToTss <- function(usefuLoci, locus2oneGene,
                                    genesColName="Gene_Symbol",
                                    coordinatesColName="Genomic_Coordinate",
                                    lociColName="probeID", doWarn=FALSE, verbose=0){
    ## Habil wrote this function to compute the distance of the input loci to the closest TSS.
    ## This function is used in plotLociTss().
    ## locus2oneGene: A data frame with the following columns:
    ## "probeID" "Gene_Symbol" "Chromosome" "Genomic_Coordinate"
    ## usefuLoci: A list named by genes. Each entry is a character vector of loci, e.g,:
    ## c("cg05116342", "cg15516558").

    result <- list()
    starTime <- Sys.time()
    message.if(paste("Computing distance to TSS for", length(usefuLoci), "genes started at:", 
                   format(starTime, usetz=TRUE), "\n"), verbose=verbose-3)
    message.if(verbose=verbose, paste("Estimated time to complete:",
                                    length(usefuLoci)/12000, "hours.\n"))
    toClosesTss <- c()
    cols1 <- c("TXSTART","TXEND","TXNAME","TXCHROM","ENTREZID") ## Not used.
    notFound <- c()
    transcriptLen <- c()
    ## require(Homo.sapiens)
    gotHomoSapiens <- get("Homo.sapiens")
    for(g1 in names(usefuLoci)){ ## Each gene,
        ind <- which(names(usefuLoci)==g1)
        if(!ind%%100)
            message.if(verbose=verbose-1, paste("Upto:  === === === > ", ind))
        ##l2g1 is a sub matrix of locus2oneGene, which is used
        ##^to identify the loci that are related to gene g1.
        l2g1 <- locus2oneGene[locus2oneGene[, genesColName]==g1, , drop=FALSE]
        l2g1 <- l2g1[!is.na(l2g1[, genesColName]),]
        l2g1 <- l2g1[l2g1[, lociColName] %in% usefuLoci[[g1]], ]
        ## Gene coordinates:
        t1 <- try(selected <- suppressMessages(AnnotationDbi::select(gotHomoSapiens, 
                                                    keys=g1, columns=cols1,
                                                    keytype="SYMBOL")), silent=verbose<2)
        if(inherits(t1, "try-error")){
            notFound <- c(notFound, g1)
            message.if(verbose=verbose-2, g1)
            next
        }
        starts <- selected[,"TXSTART"] ## A vector. length= # of transcripts of g1.
        names(starts) <- selected[,"TXNAME"]
        transcriptLen <- c(transcriptLen, selected[,"TXEND"] - starts)
        ## loci is a numeric matrix. Columns are loci, rows are the transcripts of g1, and
        ##^the values are the genomic coordinates of the loci. Repeated on columns.
        loci <- matrix(rep(l2g1[, coordinatesColName], each=length(starts)), ncol=nrow(l2g1))
        colnames(loci) <- rownames(l2g1) ## Each column corresponds to a locus.
        ## Distance to TSS:
        distanceToTss <- loci - starts ## R will repeat starts in each column.
        ## From each column, find the entry with the minimum absolute value.
        minRowInds <- Rfast::colMins(abs(distanceToTss), value=FALSE)
        ## Look at the distanceToTss marix as a vector:
        d1 <- distanceToTss[minRowInds+ nrow(distanceToTss)*(0:(ncol(distanceToTss)-1))]
        ## The length of d1 is equal to the number of transcripts of g1.
        names(d1) <- colnames(distanceToTss)
        toClosesTss <- c(toClosesTss, d1)
    }
    ##End for(g1).
    
    nas <- names(which(is.na(toClosesTss)))
    message.if(paste("Number of NAs in toClosesTss:", length(nas), "\n"),
               verbose=verbose-1)    
    result[["distanceToClosesTssPerGene"]] <- toClosesTss
    ##^ An integer vector named by the loci. If a locus is related to k genes,
    ## there would be k entries, e.g., cg27665767.4, cg27665767.8, cg27665767.9, etc.
    message.if(verbose=verbose-1, "Finding the minimum among all genes for each locus:")
    distanceToClosesTssPerlocus <- c()
    withoutDot <- gsub(names(toClosesTss), pattern="\\..*$", replacement="")
    for(l1 in unique(withoutDot)){
        min1 <- min(toClosesTss[withoutDot==l1], na.rm=TRUE)
        if(is.finite(min1)) ##If toClosesTss is NA, min1 is Inf.
            distanceToClosesTssPerlocus[l1] <- min1
    }
    result[["distanceToClosesTss"]] <- distanceToClosesTssPerlocus

    ## Output:
    message.if(verbose=verbose-1, "The transcripts length:")
    message.if(verbose=verbose-1, summary(transcriptLen))
    timeTaken <- Sys.time()-starTime
    result[["transcriptLen"]] <- transcriptLen
    result[["notFound"]] <- notFound
    if(length(notFound)>0 & doWarn)
        warning(length(notFound)," gene symbols were not found!")
    result[["timeTaken"]] <- timeTaken
    message.if(paste("Computing distance to TSS for took:", timeTaken,
                   attr(timeTaken,"units"), "\n"), verbose=verbose-3)
    return(result)
}
