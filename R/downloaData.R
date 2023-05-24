downloaData <- function(dataProject, savePath, doDnam=TRUE, doExpr=TRUE,
                        doClinical=TRUE, doMirna=FALSE, verbose=0){
    ## Download data from TCGA portal.
    ## Inputs:
    ## dataProject: A named for project name in TCGA portal.
    ## savePath: A path that data must save there.
    ## doDnam: Logical, TRUE indicates Dnam data will be downloaded.
    ## doClinical: Logical, TRUE indicates clinical data will be downloaded.
    ## doExpr: Logical, TRUE indicates expression data will be downloaded.
    ## doMirna: Logical, TRUE indicates miRNA expression data will be downloaded.
    ## ^E.g.Primary Tumor Samples
    ## Output:
    result <- list()
    result[["call"]] <- match.call()
    result[["savePath"]] <- savePath

    ## Timing:
    starTime <- Sys.time()
    message.if(me=paste("Downloading started at:", format(starTime, usetz=TRUE),"\n"), 
               verbose=verbose-2)
    message.if(me=paste("Data will be saved at:", savePath, "\n"), verbose=verbose-2)
    
    timeTaken <- list()
    queries <- list()
    timeTaken[["start"]] <- starTime
    result[["dataProject"]] <- dataProject
    ##result[["categories"]]  <- TCGAbiolinks:::getProjectSummary(project=dataProject) 
    
    dnamEndTime <- Sys.time()
    if(doDnam){
        query <- TCGAbiolinks::GDCquery(project= dataProject,
                                        data.category="DNA Methylation",
                                        platform="Illumina Human Methylation 450",
                                        data.type="Methylation Beta Value")
        ##legacy=legacy) Not supported anymore :(
        queries[["rawDnam"]] <- query
        TCGAbiolinks::GDCdownload(query, directory=savePath)
        message.if(me="Preparing DNA methylation data...", verbose=verbose+1)
        ## The sesame package is needed for the following line.
        Data <- TCGAbiolinks::GDCprepare(query, directory=savePath)
        result[["rawDnam"]] <- Data 
        dnamEndTime <- Sys.time()
        message.if(me=paste("Downloading DNA methylation data finished at:", 
                            format(dnamEndTime, usetz=TRUE), "\n"), verbose=verbose-1)
        timeTaken[["dnam"]] <- dnamEndTime-starTime
    }

    clinicalEndTime <- Sys.time()
    if(doClinical){
        query <- TCGAbiolinks::GDCquery(project=dataProject,
                                        data.category="Clinical",
                                        file.type="xml")
        queries[["clinical"]] <- query
        TCGAbiolinks::GDCdownload(query, directory=savePath)
        clinical <- TCGAbiolinks::GDCprepare_clinic(query,directory=savePath,clinical.info="patient")
        result[["clinical"]] <- clinical
        clinicalEndTime <- Sys.time()
        message.if(me=paste("Downloading clinical data finished at:",
                            format(clinicalEndTime, usetz=TRUE), "\n"), verbose=verbose-1)
        timeTaken[["clinical"]] <- clinicalEndTime-dnamEndTime
    }

    exprEndTime <- Sys.time()
    if(doExpr){
        if(FALSE){ ##  Not supported anymore :(
            query <- TCGAbiolinks::GDCquery(project=dataProject,
                                            data.category="Transcriptome Profiling",
                                            data.type="Gene expression quantification",
                                            platform="Illumina HiSeq",
                                            file.type ="normalized_results",
                                            experimental.strategy="RNA-Seq",
                                            legacy=TRUE)
        }
        query <- TCGAbiolinks::GDCquery(project=dataProject,
                                        data.category="Transcriptome Profiling",
                                        data.type="Gene Expression Quantification",
                                        workflow.type = "STAR - Counts")
        queries[["genExpr"]] <- query
        TCGAbiolinks::GDCdownload(query, directory=savePath)
        Data <- TCGAbiolinks::GDCprepare(query, directory=savePath)
        genExpr <- SummarizedExperiment::assay(Data, "fpkm_unstrand")
        offset0 <- unique(sort(genExpr))[2]
        message.if(me=paste("Log-transform after adding the smallest non-zero value:",
                            offset0), verbose=verbose-2)
        genExpr <- log10(genExpr+offset0)
        result[["genExpr"]] <- genExpr
        exprInfo <- SummarizedExperiment::colData(Data)
        result[["genExprSampleInfo"]] <- exprInfo
        exprEndTime <- Sys.time()
        message.if(me=paste("Downloading expression data finished at:",
                            format(exprEndTime, usetz=TRUE), "\n"), verbose=verbose-1)
        timeTaken[["genExpr"]] <- exprEndTime-clinicalEndTime
    }

    exprEndTime <- Sys.time()
    if(doMirna){
        query <- TCGAbiolinks::GDCquery(project=dataProject,
                                        experimental.strategy="miRNA-Seq",
                                        data.category="Transcriptome Profiling",
                                        data.type="miRNA Expression Quantification")
        queries[["mirna"]] <- query
        TCGAbiolinks::GDCdownload(query, directory=savePath)
        mirna <- TCGAbiolinks::GDCprepare(query,directory=savePath)
        result[["mirna"]] <- mirna
        mirnaEndTime <- Sys.time()
        message.if(me=paste("Downloading miRNA data finished at:",
                            format(mirnaEndTime, usetz=TRUE), "\n"), verbose=verbose-1)
        timeTaken[["mirna"]] <- mirnaEndTime-exprEndTime
    }

    timeTaken[["total"]] <- clinicalEndTime-starTime
    result[["queries"]] <- queries
    result[["timeTaken"]] <- timeTaken        

    ## Save all:
    tcga <- result
    save(tcga, file=file.path(savePath, "tcga.RData"))
    return(tcga)
}
