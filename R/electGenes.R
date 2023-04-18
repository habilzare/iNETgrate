electGenes <- function(genExpr, dnam, survival, savePath, event="Dead", locus2gene,
                        genesColName="Gene_Symbol", lociColName="probeID",
                        ratio=c(genExpr=1/3, dnam=1),
                        minCor=c(genExpr=0.2, dnam=0.2), doAddTitle=TRUE,
                        doAlLoci=FALSE, verbose=0){
    ## This is a wrapper function that executes, filterLowCor and plotLociNum 
    ## functions to filter gene expression dnam data and computes a union of
    ## selected genes.
    ## Input:
    ## dnamData: A matrix of beta values where row names are probe or loci IDs and column names 
    ## are patient IDs or case IDs.
    ## survival: A matrix where rows are patients and columns are “Time” (to event), event (i.g., "Dead")
    
    ## Start Time
    starTime <- Sys.time()
    message.if(me=paste("Cleaning all data started at:", format(starTime, usetz=TRUE)),
               verbose=verbose-3)

    result <- list()
    ## filter low cor genes
    message.if("Filtering gene expression data...", verbose=verbose)
    filteredExpr <- filterLowCor(Data=genExpr, survival=survival,
                                   whichData="expression", ratio=ratio["genExpr"],
                                   minCor=minCor["genExpr"], event=event,
                                   savePath=savePath, verbose=verbose)

    ## filter low cor loci
    message.if("Filtering DNA methylation data...", verbose=verbose)
    filteredDnam <- filterLowCor(Data=dnam, survival=survival,
                                   whichData="dnam", ratio=ratio["dnam"], 
                                   minCor=minCor["dnam"], event=event,
                                   savePath=savePath, verbose=verbose)
    
    ## Plot loci per gene:
    if(!is.null(savePath)){
        message.if("Combining filtered data...", verbose=verbose)
        plotPath <- file.path(savePath, "plots")
        dir.create(plotPath, recursive=FALSE, showWarnings=FALSE)
        plotFile <- file.path(file.path(plotPath, "pruned_loci_num_per_gene.png"))
        plotted <- plotLociNum(locus2gene=locus2gene, selectedLoci=filteredDnam$selectedList,
                                plotFile=plotFile, doAddTitle=doAddTitle)
    }

    ## Compute union:
    unionSet <- computeUnion(Genes=rownames(genExpr),
                              selectedGenes=filteredExpr$selectedList,
                              loci=rownames(dnam),
                              selectedLoci=filteredDnam$selectedList,
                              locus2gene=locus2gene, genesColName=genesColName,
                              lociColName=lociColName, doAlLoci=doAlLoci,
                              verbose=verbose)
    result[["selectedGenesExpr"]] <- filteredExpr$selectedList
    result[["selectedLoci"]] <- filteredDnam$selectedList
    result[["selectedGenesDnam"]] <- unionSet$dnamSelectGenes 
    result[["unionGenes"]] <- unionSet$unionGenes
    result[["unionLoci"]] <- unionSet$unionLoci

    ## End time:
    timeTaken <- Sys.time()-starTime
    message.if(me=paste("Filtering data took:", timeTaken, attr(timeTaken,"units"), "\n"), 
               verbose=verbose-3)

    return(result)
}
