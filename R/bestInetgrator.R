bestInetgrator <- function(bestPvalues, usefuLoci, lociPigen, netPath, verbose=0){
    ## It returns the integrator object corresponding to the best p-value from the results of the
    ##^running pipeline on the training dataset (TCGA for AML disease).
    ## Input:
    ## bestPvalues: a data frame that has three cols: "mu"=the mu value, "pvalue"=the p-value
    ##^of the related model, and "features"=the name of the model e.g., ME14.e_ME24.e_ME3.e
    ## Use a dara frame with 1 row if you want to manually choose the model.
    ## usefuLoci: A list of genes and their selected loci in the process of computing eigenloci
    ## lociPigen: A list of genes (which are the same as selectedloci’s genes) and their
    ## Related pigengene object (based on their mapped loci).
    ## netPath: Path to the network folder, where the pigengene and eigengene objects are saved
    ##^for each modules.
    ## ^For TCGA, we used networkPath[[tIndx]].
    ## besr

    ##QC:
    if(!inherits(bestPvalues, "data.frame"))
        stop("bestPvalues must be a data frame!")
    if(sum(c("features","pvalue","mu") %in% colnames(bestPvalues)) <3)
        stop("bestPvalues must have these columns: features, pvalue, and mu!")

    message.if("Making the best inetgrator...", verbose=verbose)

    ##Finding best mu based on best pvalue:
    bestPvaluesOrdered <- bestPvalues[order(bestPvalues[,"pvalue"]),]
    bestPvalue <- bestPvaluesOrdered[1,"pvalue"]
    modelName <- bestPvaluesOrdered[1,"features"]
    mu <- bestPvaluesOrdered[1,"mu"]

    message.if(me=paste("The best mu value based on p-values(",
                      bestPvalue, ") is: ", mu,"\n"), verbose=verbose-1)
    message.if(me=paste("The best modules are: ", as.character(modelName), "\n"),
               verbose=verbose-1)
    
    ## Identifying the modules for best mu
    features <- unlist(strsplit(as.character(modelName), "_"))
    moduleCombMethods <- c()
    moduleName <- c()
    pigengenes <- list() 
    for(m1 in features){
        mdl <- unlist(strsplit(m1, "\\."))
        moduleName <- c(moduleName, mdl[1]) ## e.g., "ME46"
        combMethod <- mdl[2] ## e, m or em
        moduleCombMethods <- c(moduleCombMethods, combMethod)         
        ## Getting pigengene object for best mu:
        file1 <- file.path(netPath,paste0("mu",mu),
                           paste0("Pigen_",combMethod),"pigengene.RData")
        message.if(paste("Loading pigengene data for \n", file1, "\n"),
                   verbose=verbose-2)
        pigengene <- get(load(file1, verbose=TRUE))
        pigengenes[[combMethod]] <- pigengene
    }
    ## A module-method matrix.
    moMe <- cbind(ME=moduleName, eOrM=moduleCombMethods)
    moMe <- paste(moMe[,"ME"], moMe[,"eOrM"], sep=".")
    names(moMe) <- features
  
    ## Collecting module information
    inetgrator <- computeInetgrator(mu=mu, usefuLoci=usefuLoci,
                                     genePigens=pigengenes,
                                     lociPigen=lociPigen, moMe=moMe)
    ## Plot heatmaps of the best modules: TODO
    message.if("Plotting heatmaps for selected modules...", verbose=verbose)
    Pigengene::module.heatmap(pigengene=pigengenes[[1]], c5Tree=NULL, 
                              saveDir=netPath, pos=0, 
                              mes=moduleName, verbose=verbose-1)

    message.if(paste("Finished data processing at", format(Sys.time(), usetz=TRUE),"\n"),
               verbose=verbose-3)
    return(inetgrator)
}
