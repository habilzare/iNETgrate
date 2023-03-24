sample.to.patient <- function(sampleInfo, sampleCol="barcode", patientCol="patient", 
                              sampleTypeCol="shortLetterCode",
                              sampleTypesIn=NULL, isFfpe=NULL, 
                              doRemoveDup=TRUE, verbose=0){
    ## isFfpe: A character vector of boolean values that determines in which type of samples we are interested in.
    ##^ E.g., c(FALSE, TRUE) will include all samples, TRUE will include only FFPE samples, and FALSE
    ##will exclude FFPE samples. --Habil,
    ##^ Set it to NULL to disable it.
    ##2020-09-24.
    ## Start time
    starTime <- Sys.time()
    s2p <- list()
    excludedSamples <- list()
    s2p[["doRemoveDup"]] <- doRemoveDup
    message.if(cat("Sample to patient mapping at:", format(starTime, usetz=TRUE), 
                   "\n"), verbose=verbose-3)

    ## if sampleTypeCol==NULL
    if(is.null(sampleTypeCol)){
        sample2patient <- setNames(sampleInfo[, patientCol], nm=sampleInfo[, sampleCol])
        s2p[["sample2patient"]] <- sample2patient
        return(s2p)
    }
    
    ## FFPE filtering:
    ffpeSamples <- c()
    if(!(is.null(isFfpe))){
        sampleInfo <- sampleInfo[which(sampleInfo[,"is_ffpe"] %in% isFfpe), ]
        ffpeSamples <- sampleInfo[which(sampleInfo[,"is_ffpe"] == TRUE), sampleCol]
        if(length(isFfpe) == 1 && !isFfpe){
            excludedSamples[["FFPExcluded"]] <- ffpeSamples
        }
    }
    if(nrow(sampleInfo) == 0)
        stop("Please correct isFfpe input!")

    sampleSizes <- as.matrix(table(sampleInfo[ ,sampleTypeCol]))
    sampleTypes <- rownames(sampleSizes)
    message.if(cat("Sample types in data are: \n"), verbose=verbose-1)
    message.if(print(sampleSizes), verbose=verbose-1)

    ## sampleTypes excluded
    if(is.null(sampleTypesIn)){
        message.if("All the sample types will be included in the study!",
                   verbose=verbose)
        sampleTypesIn <- sampleTypes
    } else {
        sampleTypesExcluded <- sampleTypes[!sampleTypes %in% sampleTypesIn]
        message.if(paste(sampleTypesExcluded, "sample types are excluded from
                         further analysis"), verbose=verbose)
        eSamples <- rownames(sampleInfo)[which(sampleInfo[ , sampleTypeCol] %in% sampleTypesExcluded)]
        eSamplesNamed <- setNames(sampleInfo[eSamples, sampleTypeCol], 
                                  nm=sampleInfo[eSamples, sampleCol])
        eSamplesList <- split(x=names(eSamplesNamed), f=eSamplesNamed)
        excludedSamples <- eSamplesList
    }
    
    sampleSizes <- sampleSizes[sampleTypesIn, ]
    sInfo <- sampleInfo[sampleInfo[, sampleTypeCol] %in% sampleTypesIn, 
                        c(sampleCol, patientCol, sampleTypeCol)]
    
    if(length(sampleTypesIn) > 1){
        patientSuffix <- head(names(sampleSizes)[order(sampleSizes)],
                              length(sampleSizes)-1)
        message.if(cat("Suffix to patient IDs will be added for sample types:",
                       patientSuffix, "\n"), verbose=verbose)
        sRows <- rownames(sInfo)[which(sInfo[, sampleTypeCol] %in% patientSuffix)]
        sInfo[sRows, patientCol] <- paste0(sInfo[sRows, patientCol], "_",
                                           sInfo[sRows, sampleTypeCol])
    }
   
    ## Remove any remaining duplicate samples
    if(doRemoveDup){
        dups <- find.tcga.duplicates(sampleInfo=sInfo, sampleCol=sampleCol,
                                     patientCol=patientCol,
                                     sampleTypeCol=sampleTypeCol,
                                     doTag=FALSE, verbose=verbose)
        excludedSamples[["excludeDups"]] <- dups$excludeDups
        sInfo <- sInfo[!(sInfo[, sampleCol] %in% dups$excludeDups), ]
    }

    sample2patient <- setNames(sInfo[, patientCol], nm=sInfo[, sampleCol])
    patient2type <- setNames(sInfo[, sampleTypeCol], nm=sInfo[, patientCol])

    s2p[["excludedSamples"]] <- excludedSamples
    s2p[["sample2patient"]] <- sample2patient
    s2p[["patient2type"]] <- patient2type
       
    ## End time:
    timeTaken <- Sys.time()-starTime
    message.if(cat("Sample to patient mapping took:", timeTaken, 
                   attr(timeTaken, "units"), "\n"), verbose=verbose-3)
    return(s2p)
}
