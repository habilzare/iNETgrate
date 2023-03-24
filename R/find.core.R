find.core <- function(Data, Labels,  Label1, Label2, verbose=0){
    ## Uses hierarchical clustering and a greedy approach to
    ## find "the core" (the most connected cluster of points) in the Data.
    ## Modularity criterion: average edge weight in the cluster.
    ## Data: the beta values of the loci (on columns) corresponding to 1 gene for all samples (rows).
    if(length(intersect(unique(Labels), c(Label1, Label2)))!=2)
        stop(paste("There are no patients with a label as:", Label1, "and/or",
                   Label2))
    result <- list()
    ##QC:
    checked <- check.pigengene.input(Data=Data, Labels=Labels)
    
    similarity1 <- abs(stats::cor(Data))
    g1 <- igraph::make_full_graph(n=ncol(Data))
    igraph::E(g1)$weight <- as.dist(similarity1)

    ## Only 1 locus:
    if(ncol(Data)==1){
        result[["PC1"]] <- Data##PCA on one data point is meaningless
        lociName<-colnames(Data)
        result[["lociNames"]]<-lociName
        if(verbose)
            print("n=1")
        return(result)
    }

    ## Fewer than 6 loci:
    if(ncol(Data)<6){
        communities <- rep(1,ncol(Data))
        if(verbose)
            print("n<6")
    } else {
        cfg <- igraph::cluster_fast_greedy(g1)
        communities <- igraph::membership(cfg)
    }
    weightAve <- c() ## As a measure of modularity as in FEM.
    for(c1 in unique(communities)){
        inds <- which(communities==c1)
        weightAve[c1] <- mean(similarity1[inds,inds])
    }

    ## Identify the best community:
    best <- which.max(weightAve)
    if(sum(communities==best)==1){
        ## There is only one data point in the best community
        result[["PC1"]] <- Data[,communities==best,drop=FALSE]
        lociName<-colnames(Data[,communities==best,drop=FALSE])
        result[["lociNames"]]<-lociName
        result[["aveWeight"]] <- c(weightAve)
        result[["communities"]] <- c(communities)
        result[["num"]] <- c(length(unique(communities))) 
        return(result)
    }
    
    ##modules : to identify for compute.pigengene that all loci are in the same cluster
    lociModules <- c(rep(1,sum(communities==best)))
    names(lociModules) <- colnames(Data[,communities==best,drop=FALSE])
    
    ##mainGroup: union of High or low-risk patients
    mainGroup <- Labels[c(which(Labels==Label1), which(Labels==Label2))]
    pigengene <- compute.pigengene(Data=scale(Data[names(mainGroup), communities==best,drop=FALSE]), 
                                   modules=lociModules,Labels=mainGroup,
                                   doPlot=FALSE, saveFile=NULL)
    
    projectedLoci <- project.eigen(Data=scale(Data[,communities==best,drop=FALSE]), pigengene=pigengene,
                                   naTolerance=0.05,verbose=verbose-1, ignoreModules=c())
    
    result[["PC1"]]<-as.matrix(projectedLoci$projected)
    result[["pigenObject"]]<-pigengene
    result[["lociNames"]]<-colnames(Data[,communities==best,drop=FALSE])
    result[["aveWeight"]] <- c(weightAve)
    result[["communities"]] <- c(communities)
    result[["num"]] <- c(length(unique(communities)))    
    return(result)
}

