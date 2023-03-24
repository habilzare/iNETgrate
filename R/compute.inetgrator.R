compute.inetgrator <- function(mu, usefuLoci, genePigens, lociPigen, moMe){
    ## Hanie wrote this function on 24-July-2018.
    ## It computes the integrator module.
    ## Inputs:
    ## mu: The mu value used in network construction of the training phase. It ranges 
    ##^between 0 and 1.
    ## genePigen: A list of pigengenes objects with names from c("m", "e", "em"). 
    ## usefuLoci: A list of genes and their selected loci in the process of computing eigenloci.
    ## lociPigen: A list of pigengene objects named with genes.
    ##^These genes must be a subset of the selectedloci's genes.
    ## moMe: A character vector containing selected eigengenes to be inferred
    ##^e.g. c("ME51.m", "M55.e", "M46.em")
    ## Output: 
    ## integrator: (for now a list, later an object of class "inetgrator", which we will define) 
    
    ##Check:
    if(mu<0|1<mu)
        stop("The value of mu should be between 0 and 1!")
    if(any(!names(genePigens) %in% c("m", "e", "em")))
        stop("The names of the genePigen are not defined properly!")
    if(sum(!names(lociPigen) %in% names(usefuLoci))!=0)
        stop("lociPigen's genes must be a subset of the selectedloci's genes")
    
    ##inetgrator object
    inetgrator <- list()
    inetgrator[["mu"]] <- mu
    inetgrator[["moMe"]] <- moMe
    inetgrator[["pigengenes"]] <- genePigens

    ## modules is a list of length=number of modules. For each module, 
    ## we have a weight matrix and a genelist
    modules <- list()

    ## Identifying cluster of each genes: clustering is the same for all pigengenes in the
    ##^genePigen list, so I used the first pigengene object.  
    for(n1 in seq_along((moMe))){
        moMe1 <- moMe[n1] 
        moduleName <- unlist(strsplit(moMe1, split="\\."))[1]
        eOrM1 <- unlist(strsplit(moMe1, split="\\."))[2]
        ##
        modules[[moduleName]] <- list()
        moduleNumber <- gsub(moduleName, pattern="ME", replacement="")
        geneColor <- genePigens[[1]]$orderedModules
        moduleGenes <- names(geneColor[which(geneColor==moduleNumber)])
        geneWeights <- genePigens[[eOrM1]]$weights[moduleGenes, "Weight"]
        modules[[n1]][["weights"]] <- cbind(modules[[n1]][["weights"]], geneWeights)
        colnames(modules[[n1]][["weights"]])[ncol(modules[[n1]][["weights"]])] <- eOrM1
        genes <- list()
        for(g1 in moduleGenes){
            genes[[g1]] <- list()
            if(length(usefuLoci[[g1]])==0)
                stop(paste("No locus was founded for gene",g1,"!"))
            if(length(usefuLoci[[g1]])>1){
                genes[[g1]][["loci"]] <- as.matrix(lociPigen[[g1]]$membership)
                genes[[g1]][["pigengene"]] <- lociPigen[[g1]]
            } else {##lenght is 1
                genes[[g1]][["loci"]] <- matrix(data=1,dimnames=list(usefuLoci[[g1]]))
            }
                colnames(genes[[g1]]$loci) <- "membership"
            } ## for(g1 in moduleGenes)
        modules[[n1]][["genes"]] <- genes ##contains eigenloci info
    }##for(n1 in names(moduleInfo))
    
    inetgrator[["modules"]] <- modules

    ## return
    return(inetgrator)   
}    
