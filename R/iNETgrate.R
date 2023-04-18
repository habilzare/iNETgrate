iNETgrate <- function(Data, clinSettings, mus=(0:10)/10, saveDir="iNETgrate", 
                      annLib="Auto", isFfpe=c(FALSE), doRemoveTOM=TRUE, 
                      minModuleSize=5, corMethod="pearson", doReturNetworks=FALSE, 
                      RsquaredCut=0.75, combiningMu=NA, favRisk="High", 
                      subSet="Int", xmax1=15, xmin1=0, doCox=TRUE, eOrMs="e", 
                      time2day=1, until=1, minRecall4L=0.2, 
                      minRecall4H=0.05, verbose=0){

    ## idType="ENTREZID", pathwayDb=NULL, OrgDb=org.Hs.eg.db){
   
    results <- list()
    results[["call"]] <- match.call()
    Label1 <- "Low"
    Label2 <- "High" 

    ## QC
    if(!inherits(Data, "list"))
        stop("Data must be a named list!")
    if(length(Data) < 4 & length(names(Data)) < 4)
        stop("Data should have minimum 4 components and their names!")
    expecteData <- c("genExpr", "genExprSampleInfo", "rawDnam", "clinical")
    if(any(! expecteData %in% names(Data)))
        stop("Please refer to documentation for expected component names of Data!")
    if(!inherits(clinSettings, "character"))
        stop("clinSettings must be a named vector!")
    cliNames <- c("patientIDCol","eventCol", "timeCol", "riskCatCol",
                  "riskFactorCol", "event", "otherLabel", "riskHigh", "riskLow")
    if(any(! cliNames %in% names(clinSettings)))
        stop("Please refer to documentation for expected names of clinSettings!")
    ##if(length(clinSettings) != length(cliNames))
    ##    stop("Please add names for every component of the vector!")
    ##^ More names is OK, Habil.
    otherLabel <- clinSettings["otherLabel"]
    if(otherLabel=="NULL")
        otherLabel <- NULL
    if(length(grep(saveDir, pattern=" ")>0))
        stop("saveDir cannot have space!")

    ## Save results
    dir.create(path=saveDir, recursive=TRUE, showWarnings=FALSE)

    ## Clean data:
    message.if("Step 1/7: Cleaning all the Data...", verbose=verbose)
    cleaned <- cleanAllData(genExpr=Data$genExpr, 
                             genExprSampleInfo=Data$genExprSampleInfo, 
                             rawDnam=Data$rawDnam, savePath=saveDir, 
                             annLib=annLib, clinical=Data$clinical, 
                             riskCatCol=clinSettings["riskCatCol"], 
                             riskFactorCol=clinSettings["riskFactorCol"], 
                             otherLabel=otherLabel, 
                             riskHigh=clinSettings["riskHigh"], 
                             riskLow=clinSettings["riskLow"], 
                             patientIDCol=clinSettings["patientIDCol"], 
                             eventCol=clinSettings["eventCol"],
                             event=clinSettings["event"], 
                             timeCol=clinSettings["timeCol"], verbose=verbose-1)

    ## Select genes
    message.if("Step 2/7: Filtering Data...", verbose=verbose)
    elected <- electGenes(genExpr=cleaned$genExpr, dnam=cleaned$dnam,
                           survival=cleaned$survival, savePath=saveDir, 
                           locus2gene=cleaned$locus2gene, 
                           doAlLoci=FALSE, verbose=verbose-1)
    results[["unionGenes"]] <- elected$unionGenes

    ## Computing eigenloci
    message.if("Step 3/7: Computing eigenloci...", verbose=verbose)
    patientLabel <- setNames(as.character(cleaned$survival[,"Risk1"]),
                             nm=rownames(cleaned$survival))
    inBoth <- intersect(colnames(cleaned$dnam), names(patientLabel))
    LabelsIn <- patientLabel[names(patientLabel) %in% inBoth]

    message.if(paste("Dropping ", length(patientLabel)-length(inBoth), "patients\n",
                   "because of missing survival, expression or methylation data\n"), 
               verbose=verbose)
    results[["Labels"]] <- LabelsIn

    computedEloci <- computeEigenloci(dnam=cleaned$dnam[ ,inBoth], 
                                       geNames=elected$unionGenes,
                                       locus2gene=cleaned$locus2gene, 
                                       Labels=LabelsIn, plotPath=saveDir,
                                       Label1=Label1, Label2=Label2,
                                       dnamGene=NULL, doDebug=FALSE, 
                                       verbose=verbose-1)

    eigenloci <- computedEloci$eigenloci
    results[["eigenloci"]] <- eigenloci

    ## Make combined network
    message.if("Step 4/7: Making a integrative network...", verbose=verbose)
    netPath <- file.path(saveDir, "network")
    dir.create(netPath, showWarnings=FALSE, recursive=TRUE)
    message.if(paste("Network data is saved at: \n", netPath, "\n"), verbose=verbose)

    madeNetwork <- makeNetwork(genExpr=cleaned$genExpr, eigenloci=eigenloci,
                                geNames=elected$unionGenes, mus=mus, 
                                doRemoveTOM=doRemoveTOM, outPath=netPath, 
                                minModuleSize=minModuleSize, corMethod=corMethod,
                                doReturNetworks=doReturNetworks, 
                                RsquaredCut=RsquaredCut, verbose=verbose-1)

    ## Get Eigengenes
    message.if("Step 5/7: Computing eigengenes...", verbose=verbose)
    eigenGenes <- computeEigengenes(genExpr=cleaned$genExpr, eigenloci=eigenloci, 
                                      netPath=netPath, geNames=elected$unionGenes,
                                      Labels=patientLabel, Label1=Label1, 
                                      Label2=Label2, mus=mus,
                                      combiningMu=combiningMu, 
                                      survival=cleaned$survival, 
                                      event=clinSettings["event"], 
                                      verbose=verbose-1, 
                                      mu2modules=madeNetwork$mu2modules)

    ## Survial Analysis
    message.if("Step 6/7: Survial analysis for test set...", verbose=verbose)
    survivalPath <- file.path(netPath, "survival")
    dir.create(survivalPath, showWarnings=FALSE, recursive=TRUE)
    message.if(paste("Survival output is saved at: \n", survivalPath, "\n"), 
               verbose=verbose)

    survivalAnalysed <- analyzeSurvival(survival=cleaned$survival, 
                                          Labels=patientLabel, favRisk=favRisk, 
                                          subSet=subSet, mus=mus, 
                                          netPath=netPath, outPath=survivalPath,
                                          doCox=doCox, eOrMs=eOrMs,
                                          time2day=time2day, until=until,
                                          xmax1=xmax1, xmin1=xmin1,
                                          minRecall4L=minRecall4L, 
                                          minRecall4H=minRecall4H, 
                                          verbose=verbose-1)

    ## Indentify best iNETgrator
    message.if("Step 7/7: Identifying the best module(iNETgrator)...", verbose=verbose)
    message.if(paste("Best inetgrator output is saved at: \n", netPath, "\n"), 
               verbose=verbose)

    bestPval <- read.csv(file=file.path(survivalPath, "bestPvalues_e.csv"), header=TRUE)
    bestPval <- bestPval[, -1]
    inetgrator <- bestInetgrator(bestPvalues=bestPval, usefuLoci=computedEloci$usefuLoci, 
                                  lociPigen=computedEloci$lociPigen, netPath=netPath)
    results[["bestInet"]] <- inetgrator

    return(results)
}
