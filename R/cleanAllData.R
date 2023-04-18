cleanAllData <- function(genExpr, genExprSampleInfo, rawDnam, rawDnamSampleInfo=NULL, 
                          sampleCol="barcode", patientCol="patient", 
                          sampleTypeCol="shortLetterCode", sampleTypesIn=NULL,  
                          savePath, isFfpe=c(FALSE), annLib="Auto", clinical, 
                          patientIDCol="bcr_patient_barcode", eventCol="vital_status", 
                          event="Dead", timeCol="days_to_last_followup", riskFactorCol, 
                          riskCatCol=NULL, riskHigh=NULL, riskLow=NULL,
                          otherLabel=NULL, verbose=0){
    ## This is a wrapper function that executes preprocessDnam,
    ## prepareSurvival functions, and sets sample2patient mapping for geneExpr,
    ## rawDnam, and clinical data.
    ## Input:
    ## genExpr: A gene expression matrix where rows are genes and columns are
    ## samples 
    ## genExprSampleInfo: A dataframe of gene expression sample information
    ## where rownames are sample IDs and has minimum three columns with sample ID,
    ## patient ID, and sample type.
    ## sampleCol: Name of column of sample information dataframe that contains
    ## sample IDs.
    ## patientCol: Name of column of sample information dataframe that contains
    ## patient IDs.
    ## sampleTypeCol: Name of column of sample information dataframe that contains
    ## sample type information.
    ## rawDnam: Summarized experiment object with DNA methylation data similar 
    ## to TCGABiolinks DNA methylation output.
    ## rawDnamSampleInfo: A dataframe of rawDnam sample information
    ## clinical: A dataframe of clinical data where rownames are patient IDs,
    ## and containing survival data and riskfactor information at the least.
    ## riskCatCol: The column in the clinical that has the risk levels.
    ## riskFactorCol: The column in the clinical that has description of the risk (e.g., cytogenetic
    ## abnormalities in AML).
    ## otherLabel: a vector with names the same as row names of clinical.
    ## Output
    ## A list containing
    ## survival: A matrix where rows are patients and columns are “Time” (to event), event (i.g., "Dead")

    ## Start Time
    starTime <- Sys.time()
    message.if(me=paste("Cleaning all data started at:", format(starTime, usetz=TRUE), "\n"),
               verbose=verbose-3)
    message.if(paste("Cleaned data will be saved at:", 
                   savePath, "\n"), verbose=verbose-1)
    
    result <- list()
    sample2patient <- list()
    patient2type <- list()

    ## preprocess dnam
    message.if(me="Cleaning and processing dnam data...", verbose=verbose)
    processedDnam <- preprocessDnam(rawDnam=rawDnam,
                                     rawDnamSampleInfo=rawDnamSampleInfo, 
                                     savePath=savePath,  annLib=annLib,
                                     verbose=verbose)

    ## prepare survival data
    message.if(me="Processing clinical data for analyzeSurvival...",
               verbose=verbose)
    survival <- prepareSurvival(clinical=clinical, patientIDCol=patientIDCol,
                                 eventCol=eventCol, event=event,
                                 timeCol=timeCol, riskCatCol=riskCatCol, 
                                 riskFactorCol=riskFactorCol, riskHigh=riskHigh,
                                 riskLow=riskLow, otherLabel=otherLabel,
                                 verbose=verbose)

    ## sample to patient mapping
    message.if(me="Sample to patient mapping for gene expression data ...",
               verbose=verbose)
    s2pExpr <- sample2atient(sampleInfo=genExprSampleInfo, sampleCol=sampleCol,
                                 patientCol=patientCol, sampleTypeCol=sampleTypeCol,
                                 sampleTypesIn=sampleTypesIn, isFfpe=isFfpe, 
                                 verbose=verbose)
    sample2patient[["genExpr"]] <- s2pExpr$sample2patient
    if(!is.null(sampleTypeCol))
        patient2type[["genExpr"]] <- s2pExpr$patient2type
    ##^Samples that are not in sampleTypesIn are excluded.
    ## Subset genExpr data based on included samples.
    genExpr <- genExpr[ ,colnames(genExpr) %in% names(s2pExpr$sample2patient)]

    ## Changing column names to patient IDs with sample type tags
    colnames(genExpr)  <- s2pExpr$sample2patient[match(colnames(genExpr),
                                                 names(s2pExpr$sample2patient))]

    message.if(me="Sample to patient mapping for dna methylation data ...",
               verbose=verbose)
    betaMatrix <- processedDnam$dnam
    s2pDnam <- sample2atient(sampleInfo=processedDnam$sampleInfo, 
                                 sampleCol=sampleCol, patientCol=patientCol, 
                                 sampleTypeCol=sampleTypeCol, isFfpe=isFfpe,
                                 sampleTypesIn=sampleTypesIn, verbose=verbose)
    sample2patient[["dnam"]] <- s2pDnam$sample2patient
    if(!is.null(sampleTypeCol))
        patient2type[["dnam"]] <- s2pDnam$patient2type
    ##^Samples that are not in sampleTypesIn are excluded.
    ## Subset betamatrix data based on included samples.
    betaMatrix <- betaMatrix[ ,colnames(betaMatrix) %in% names(s2pDnam$sample2patient)]
    ## Changing column names to patient IDs with sample type tags
    colnames(betaMatrix)  <- s2pDnam$sample2patient[match(colnames(betaMatrix),
                                                    names(s2pDnam$sample2patient))]

    result[["survival"]] <- survival
    result[["genExpr"]] <- genExpr
    result[["sample2patient"]] <- sample2patient
    result[["patient2type"]] <- patient2type
    result[["dnam"]] <- betaMatrix
    result[["locus2gene"]] <- processedDnam$locus2gene

    ## End time:
    timeTaken <- Sys.time()-starTime
    message.if(me=paste("Cleaning all data took:", timeTaken, attr(timeTaken,"units"), "\n"), 
               verbose=verbose-3)

    ## Save data
    cleanedData <- result
    save(cleanedData, file=file.path(savePath, "cleanData.RData"))
    return(result)
}
