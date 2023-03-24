find.tcga.duplicates <- function(sampleInfo, sampleCol="barcode", 
                                     patientCol="patient", 
                                     sampleTypeCol="shortLetterCode", 
                                     doTag=FALSE, verbose=0){
    sInfoDups <- sampleInfo

    ## QC:
    if(is.null(sampleTypeCol) & doTag==TRUE)
        stop("sampleTypeCol cannot be NULL if doTag is TRUE!")
    ## doTag adds sample type tag to the patient ids. Thus avoiding removing any
    ## samples coming from different sample types.
    if(doTag){
        message.if("Patient column are tagged with sample type", verbose=verbose-2)
        sInfoDups[, patientCol] <- paste0(sInfoDups[, patientCol], "_",
                                          sInfoDups[, sampleTypeCol])     
    }

    ## Identify duplicate patients
    duPatients <- sInfoDups[which(duplicated(sInfoDups[, patientCol])),
                                patientCol]
    
    ## Identify duplicate samples
    dupSamples <- sInfoDups[sInfoDups[,patientCol] %in% duPatients, 
                                sampleCol]
    
    dupInfo <- setNames(sInfoDups[dupSamples,sampleCol],
                        nm=sInfoDups[dupSamples, patientCol])

    message.if("Duplicate samples identified are:", verbose=verbose-3)
    message.if(dupInfo, verbose=verbose-3)

    keptDups <- c()
    excludeDups <- c()
    for(p1 in duPatients){
        samples <- dupInfo[names(dupInfo)==p1]
        sorted <- sort(samples, decreasing=TRUE) ##lexicographical sorting
        keptDups <- c(keptDups, sorted[1])
        excludeDups <- c(excludeDups, sorted[-1])
    }

    message.if("Duplicate samples excluded are:", verbose=verbose-3)
    message.if(excludeDups, verbose=verbose-3)

    sampleInfo <- sampleInfo[!sampleInfo[, sampleCol] %in% excludeDups, ]
    return(list(sampleInfo=sampleInfo, keptDups=keptDups, excludeDups=excludeDups))
} 
