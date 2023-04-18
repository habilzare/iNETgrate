computeUnion <- function(Genes, selectedGenes, loci, selectedLoci, locus2gene,
                          doAlLoci=FALSE, genesColName="Gene_Symbol", 
                          lociColName="probeID", verbose=0){
    ## Compute the list of genes and loci that are correlating with survival time based on
    ##^expression levels or Dnam data

    ##Timing:
    starTime <- Sys.time()
    message.if(me=paste("Identifying a union gene set:", format(starTime, usetz=TRUE)),
               verbose=verbose-3)

    ## QC:
    if(length(intersect(locus2gene[ ,genesColName], selectedGenes))==0){
        stop("selectedGenes do not overlap with ", genesColName,
             " column in locus2gene! Are gene IDs similar?")
    }
    if(is.null(rownames(locus2gene)))
        stop("locus2gene must have row names!")  
    if(any(rownames(locus2gene) != locus2gene[, lociColName]))
        stop("Row names of locus2gene do not agree with its ", lociColName," column!")

    ## One locus - one gene mapping
    createdLocusGene <- createLocusGene(locus2gene=locus2gene,
                                         genesColName=genesColName,
                                         lociColName=lociColName, verbose=0) 
    locus2oneGene <- createdLocusGene$locus2oneGene

    ## union of gene names that have expr or dnam good data.
    lociInBeta <- rownames(locus2oneGene)[locus2oneGene[, lociColName] %in% loci]
    ##^ This line gets all the loci in preprocessed dnam data
    dnamGene <- locus2oneGene[lociInBeta, genesColName]
    dnamGene <- unique(dnamGene[which(dnamGene != "")])
    ##^ dnamGene is a list of genes that is associated with preprocessed, but
    ##NOT filtered dnam based on high correlation with survival time and vital status.

    ## Computing union for selected loci and genes
    selectedLociRows <- rownames(locus2oneGene)[locus2oneGene[, lociColName] %in%
                                               selectedLoci]
    dnamSelectGene <- locus2oneGene[selectedLociRows, genesColName]
    dnamSelectGene <- unique(dnamSelectGene[which(dnamSelectGene != "")])
     ##^ The list of genes that have at least one locus in selectedLoci.
     ## We add them to the list of already selected genes:
    unionGenes <- union(selectedGenes, dnamSelectGene)
    if(doAlLoci){
        message.if(me="All loci included in computing unionGenes",
                   verbose=verbose)
        unionGenes <- union(selectedGenes, dnamGene)
    }
    ## Exclude the genes for which we do not have expr data:
    unionGenes <- intersect(unionGenes, Genes) 
    ## Exclude the genes for which we do not have dnam data:
    unionGenes <- intersect(unionGenes, locus2oneGene[lociInBeta, genesColName])

    message.if(paste("unionGenes set size:", length(unionGenes), "\n"),
               verbose=verbose)
    
    ##Identify unionLoci: the loci that correspond to any of unionGenes
    ##AND they have valid dnam data:
    inds <- which(locus2oneGene[, genesColName] %in% unionGenes)
    unionLoci <- unique(as.character(locus2oneGene[inds, lociColName]))
    unionLoci <- unionLoci[unionLoci %in% loci]
    if(length(unionLoci)==0)
        stop("No loci left!")

    # End time:
    timeTaken <- Sys.time()-starTime
    message.if(me=paste("Computing union set took:", timeTaken, 
                      attr(timeTaken,"units"), "\n"), verbose=verbose-3)
    gc()

    ## Output:
    return(list(unionLoci=unionLoci, 
                unionGenes=unionGenes, 
                dnamSelectGenes=dnamSelectGene,
                timeTaken=timeTaken))
}
