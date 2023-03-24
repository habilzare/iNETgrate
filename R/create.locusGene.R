create.locusGene <- function(locus2gene="Auto", genesColName="Gene_Symbol", 
                             lociColName="probeID", verbose=0){
    ## Creates a dataframe useful for mapping probe IDs to genes (one to one) and vice versa.
    ## Input: preprocessed locus2gene dataframe.
    ## Output: A list including a locusGene matrix with at least "probeID" (e.g., cg00050873)
    ## and "Gene_Symbol" (gene symbol) columns.
    ## Tokenizing: 
    ## locus2gene[,"Gene_Symbol"] contains ; as seprators between gene names corresponding to the same locus.
    ## The following lines will create  a new matrix, locus2oneGene, which has seprated gene symbols.
    ## The columns are the same as locus2gene.
    ## The rownames are modified "probeID" (e.g. cg00001) so that they are unique.
    ##E.g., there are 2 rows for cg00002930:
    ##                       probeID      Gene_Symbol   Chromosome      Genomic_Coordinate
    ##cg00002930             cg00002930     NFKBIL1          6           31515398
    ##cg00002930.1           cg00002930    ATP6V1G2          6           31515398

    result <- list()
    if("Auto" %in% locus2gene){
        annLib <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
        ann450k <- minfi::getAnnotation(annLib)
        locus2gene <- ann450k[ ,c("chr", "pos", "Name", "UCSC_RefGene_Name")]
        colnames(locus2gene) <- c("Chromosome", "Genomic_Coordinate", "probeID","Gene_Symbol")
    }

    locus2gene[is.na(locus2gene[, genesColName]), genesColName] <- ""
    message.if(cat("Analyzing", nrow(locus2gene),
                    "loci, each possibly corresponding to multiple genes...\n"),
               verbose=verbose)

    ## Consider the column with gene names:
    class(locus2gene[, genesColName]) <- "character"

    ## Using tidyr.
    locus2oneGene <- tidyr::separate_rows(data=locus2gene, sep=";",
                                          tidyselect::all_of(genesColName))
    locus2oneGene <- as.data.frame(locus2oneGene)
    locus2oneGene <- unique(locus2oneGene)
    rownames(locus2oneGene) <- make.names(names=locus2oneGene[, lociColName], 
                                          unique=TRUE)

    message.if(cat("Now we have", nrow(locus2oneGene), "loci-gene pairs. \n"),
               verbose=verbose)

    result[["noGenes"]] <- rownames(locus2oneGene)[which(locus2oneGene[,genesColName]=="")]
    result[["locus2oneGene"]]  <- locus2oneGene
    return(result)
}
