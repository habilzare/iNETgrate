filter.low.cor <- function(Data, survival, whichData=c("expression", "dnam"), 
                           ratio=1/3, minCor=0.2, event="Dead", savePath,
                           verbose=0){
    ## Computes the correlation between gene expression/DNAm and survival time,
    ## and vital status and then excludes those genes that have low
    ## correlation with one/both.

    starTime <- Sys.time()
    message.if(me=cat("Cleaning all data started at:", format(starTime, usetz=TRUE)),
               verbose=verbose)
    message.if(me=cat("Dim of input Data to filter.low.cor():", dim(Data), "\n"), 
               verbose=verbose-1)
    message.if(me=cat("Correlation threshold: +/-", minCor, "\n"),
               verbose=verbose-1)
    message.if(me=cat("Keep", ratio * 100, "% of top correlated Data. \n"),
               verbose=verbose-1)

    ## Check input values:
    if(minCor==0 & ratio==1)    
        warning("All cleaned Data matrix will be returned since no threshold is defined!!")
    if(minCor!=0 & ratio!=1)
        message.if(paste("Both minCor and ratio are passed as inputs,",
                         "both thresholds will be used."), verbose=verbose)
    if(length(whichData) !=1)
        stop("Please chose correct input data type!!(expression or dnam)")
    if(minCor>1||minCor<0)
        stop("The minCor value should be between 0 and 1!!")
    if(ratio>1||ratio<0)
        stop("The ratio value should be between 0 and 1!!")

    filtered <- list()
    alList <- rownames(Data)
    ## Mapping sample to event (Dead)
    allSamples <- colnames(Data)
    deadPatients <- rownames(survival[(survival[, event]==1),])
    deadSamples <- allSamples[gsub("_.*","",allSamples) %in% deadPatients]
    deathTime <- as.numeric(survival[gsub("_.*","",deadSamples), "Time"])
    names(deathTime) <- deadSamples

    ## keep dead patients and their death time (Removing NAs):
    deathTime <- deathTime[!is.na(deathTime)]
    
    ## cor between gene expression and survival time
    common <- intersect(colnames(Data), names(deathTime))

    ## Check if patient ids and sample ids are correctly mapped:
    if(length(common)==0){
        stop(paste("Columns of",whichData,"data are different from rows of survival for", 
                   event, "patients!"))
    }       
    corTime <- stats::cor(t(Data[, common]), as.numeric(deathTime[common]),
                          method="spearman", use="everything")
    colnames(corTime) <- "Time"
    naCorTime <- sum(is.na(corTime))
    message.if(cat(naCorTime, "entities with NA correlation to survival \n", 
                   "time are removed. \n"), verbose=verbose)
    corTime <- na.omit(corTime)

    ##minCor threshold
    selectedTime <- rownames(corTime)[which(abs(corTime)>=minCor)]
    message.if(cat(length(selectedTime), "entities in", whichData, "are selected \n",
                   "as their absolute correlation with survival time >", minCor,
                   "\n"), verbose=verbose)

    ## Correlated with vital status (Dead)
    sample2vital <- setNames(survival[gsub("_.*","", allSamples), event], 
                             nm=allSamples)

    if(any(is.na(sample2vital)))
        sample2vital <- sample2vital[!is.na(sample2vital)]
    message.if(cat(length(allSamples)-length(sample2vital), " cases are \n", 
                   "removed because of NA in vital status \n"), verbose=verbose)

    ## The correlation between loci and vital status:
    corVital <- stats::cor(t(Data[, names(sample2vital)]), sample2vital,
                           method="spearman", use="everything")
    colnames(corVital) <- "Vital"
    naCorVital <- sum(is.na(corVital))
    message.if(cat(naCorVital, "entities with NA correlation to vital status \n",
                   "are removed. \n"), verbose=verbose)
    corVital <- na.omit(corVital)

    ## Identify loci that are highly correlated with vital status:
    selectedVital <- rownames(corVital)[which(abs(corVital) >= minCor)]
    message.if(cat(length(selectedVital), "entities in", whichData, "are \n",
                   "selected as their absolute correlation with vital status >", 
                   minCor, "\n"), verbose=verbose)

    ## Chose genes/loci based on either corTime or corVital or both
    corsDF <- merge(as.data.frame(corTime), as.data.frame(corVital), by=0, all=TRUE)
    cors <- as.matrix(corsDF[,c("Time", "Vital")])
    rownames(cors) <- corsDF[, "Row.names"]

    selectedCor <- union(selectedTime, selectedVital)
    message.if(cat("Number of selected data in", whichData, "based on \n", 
                   "correlation threshold", minCor, "for vital status AND \n",
                   "survival time:", length(selectedCor), ". \n"), verbose=verbose)

    ##ratio threshold:
    idx <- length(selectedCor) * ratio
    selectedCor <- selectedCor[seq_len(idx)]
    message.if(cat("Keeping only", length(selectedCor), "entities in", whichData, 
                   "because ratio is set to", ratio, ". \n"), verbose=verbose)

    if(!is.null(savePath)){
        plotPath <- file.path(savePath, "plots")
        dir.create(plotPath, showWarnings=FALSE)
        message.if(cat("The correlation plots are saved at: \n", plotPath, "\n"),
                   verbose=verbose)
        plotCors(corTimeData=corTime, corVitalData=corVital, corUnion=cors,
                 savePath=plotPath, corCutoff=minCor, whichData=whichData)
    }

    ## return results
    filtered[["alList"]] <- alList
    filtered[["selectedList"]] <- selectedCor
    filtered[["corValues"]] <- cors

    ## End time:
    timeTaken <- Sys.time()-starTime
    message.if(me=cat("Filtering data took:", timeTaken, attr(timeTaken,"units"), "\n"), 
               verbose=verbose-3)
    
    return(filtered)
}
