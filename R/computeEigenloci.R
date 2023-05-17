computeEigenloci <- function(dnam, locus2gene, geNames=NULL, Labels,
                              Label1, Label2, plotPath=NULL, genesColName="Gene_Symbol",
                              coordinatesColName="Genomic_Coordinate",
                              lociColName="probeID",dnamGene=NULL,
                              doDebug=FALSE, verbose=0){
    ## compute the eigenloci (i.e., eigenvector, PC1) of the most correlated loci.
    ## Inputs:
    ## dnam: A matrix of beta values with loci on rows and cases on columns.
    ## geNames: A character vector of selected genes that have correlation with
    ##survival time based on Expression levels or Dnam data.
    ## locus2oneGene: A matrix as explained in computeUnion.R.
    ## Labels: A character vector that determines the risk group for each patient.
    ##^It's needed to balance samples when computing the eigenloci.
    ## plotPath: A string that point to location of plots.
    ## Output:
    ## eigenloci: A matrix with samples on rows, and genes on columns in the same order as "genes"
    ## in the inputs.
    ## usefuLoci: A list of vectors.
    ##the corresponding selected loci.
    ## lociPigen: A list of Pigengene objects. Names are genes.
    ## distanceToTss: A list of distance of loci to Tss.
    ## distanceToTssDnam: A list of distance of locai to TssDnam.

    result <- list()
    ##QC:
    if(is.null(plotPath))
        stop("plotPath is null!")
    for(l1 in c(Label1, Label2)){
        if(! l1 %in% Labels)
            stop("There is no patient labeled as: ", l1)
    }
    checked <- check.pigengene.input(Data=t(dnam), Labels=Labels)

    if(is.null(rownames(locus2gene)))
        stop("locus2gene must have row names!")  
    if(any(rownames(locus2gene) != locus2gene[, lociColName]))
        stop("Row names of locus2gene do not agree with its ",
                   lociColName," column!")
         
    ## One locus - one gene mapping
    createdLocusGene <- createLocusGene(locus2gene=locus2gene,
                                         genesColName=genesColName,
                                         lociColName=lociColName, verbose=0) 
    locus2oneGene <- createdLocusGene$locus2oneGene

    if(is.null(geNames))
        geNames <- unique(locus2oneGene[rownames(dnam), genesColName])

    message.if(me="Computing gene2locus...", verbose=verbose)
    gene2loci <- list()
    for(g0 in geNames){
        gene2loci[[g0]] <- locus2oneGene[which(locus2oneGene[, genesColName] == g0), lociColName]
        gene2loci[[g0]] <- intersect(gene2loci[[g0]], rownames(dnam))
    }
 
    numG <- length(geNames)
    message.if(paste("For every one of the", numG, "genes... \n"), verbose=verbose)
    eigenloci <- c()
    lociPigen <- list()
    usefuLoci <- list()
    ##compute the eigenloci (weighted average over some selected loci) for each gene.
    for(i1 in seq_len(numG)){
        if((i1 %% 100) == 0)
            message.if(paste(i1, "of", numG, "\n"), verbose=verbose)
        g1 <- geNames[i1]
        loci1 <- gene2loci[[g1]]
        if(length(loci1)==0)
            stop("No locus is available for ", g1)
        if((i1 %% 100) ==0){
            message.if(paste(g1, i1, "th out of ", numG, ". Number of its loci:", 
                           length(loci1), "\n"), verbose=verbose)
        }
        ##trim
        found <- findCore(Data=t(dnam[loci1, , drop=FALSE]), Labels=Labels, 
                           Label1=Label1, Label2=Label2)
        eigenloci <- cbind(eigenloci, found[["PC1"]])
        colnames(eigenloci)[ncol(eigenloci)] <- g1
        usefuLoci[[g1]] <- found[["lociNames"]]        
        if(length(found[["lociNames"]]) != 1)
            lociPigen[[g1]] <- found[["pigenObject"]]
    }
    if(doDebug){
        save(eigenloci, lociPigen, usefuLoci,
             file=file.path(plotPath, "trim_debug.RData"))
        ##browser()
    } 

    ## Plot the loci position with respect to TSS:
    ## fileDistanceToTss <- gsub(plotPath, pattern="\\.RData", replacement="_distanceToTss.RData")
    fileDistanceToTss <- file.path(plotPath,"usefuLoci_distanceToTss.RData")
    plotFileGene <- file.path(plotPath,"usefuLoci_genes.png")
    distanceToTss <- distanceToTss(usefuLoci=usefuLoci,
                                             locus2oneGene=locus2oneGene,
                                             genesColName=genesColName,
                                             lociColName=lociColName,
                                             coordinatesColName=coordinatesColName,
                                             verbose=verbose-1)
    plotLociTss(distanceToTss=distanceToTss$distanceToClosesTss, plotFile=plotFileGene)

    ## Result:
    result[["eigenloci"]] <- eigenloci
    result[["usefuLoci"]] <- usefuLoci
    result[["lociPigen"]] <- lociPigen
    result[["distanceToTss"]] <- distanceToTss
    
    ## Only the genes with survival-correlating dnam, i.e., excluding the genes with
    ##survival-correlating expr:
    if(!is.null(dnamGene)){
        fileDistanceToTssDnam <- file.path(plotPath,"usefuLoci_distanceToTss_dnam.Rdata")
        usefuLociDnam <- usefuLoci[intersect(names(usefuLoci), dnamGene)]
        distanceToTssDnam <- distanceToTss(usefuLoci=usefuLociDnam,
                                                     locus2oneGene=locus2oneGene,
                                                     genesColName=genesColName,
                                                     lociColName=lociColName,
                                                     coordinatesColName=coordinatesColName)
        plotFileDnam <- file.path(plotPath,"usefuLoci_dnamGenes.png")
        plotLociTss(distanceToTss=distanceToTssDnam$distanceToClosesTss, plotFile=plotFileDnam)
        result[["distanceToTssDnam"]] <- distanceToTssDnam
    }

    ## Output:
    return(result)
}
