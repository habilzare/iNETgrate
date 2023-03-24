make.network <- function(genExpr, eigenloci, geNames, mus, doRemoveTOM=TRUE,
                         outPath, minModuleSize=5, corMethod="pearson",
                         doReturNetworks=FALSE, RsquaredCut=0.75,
                         doSaveCombined=FALSE, verbose=0){
    ##  Use adjacency to create network from input data. Then,
    ##Use mu value for create combine network from exprNetwork and dnamNetwork.
    ## Inputs:
    ## genExpr: A matix with genes on rows, and patient (tagged if needed) on columns.
    ## eigenloci: A matrix with samples on rows, and genes on columns.
    ## geNames: A character vector of selected genes.
    ## mus: A vector of mu values.
    ## doRemoveTOM: A boolean that show the TOM file must remove or not.
    ## outPath: A string to showed the path for save plot and combinednetwork.
    ## minModuleSize: The value that controls the minimum number of genes per module.
    ## corMethod: "spearman" or "pearson".
    ## doReturNetworks: A boolean value to determine whether returning exprNetwork and dnamNetwork.
    ## Output:
    ## powerVector: A numeric vector of power values named with the corresponding mus.
    ##^See WGCNA::pickSoftThreshold documentation.
    ## outliersNumber:  A named numeric vector of number of outliers per mu value.
    ## exprNetwork: A matrix that showed adjacency between gene expresion.
    ## dnamNetwork: A matrix that showed adjacency between gene methylation.
    ## library(WGCNA)

    result <- list()
    ##QC:
    if(is.null(geNames))
        stop("geNames cannot be NULL!")

    ## if(any(rownames(genExpr)!=colnames(eigenloci))){
    genExpr <- genExpr[geNames,]
    eigenloci <- eigenloci[, geNames]
        ## }
 
    if(corMethod == "pearson"){
        corInput<-paste0("use='p', method='", corMethod, "'")
    } else{
        corInput<-paste0("method='", corMethod, "'")
    }

    ## message.if("Making the integrative network....", verbose=verbose)
    message.if(paste("genExpr is of dim:", paste(dim(genExpr), collapse="*")),
               verbose=verbose-2)    
    message.if(paste("eigenloci is of dim:", paste(dim(eigenloci), collapse="*")),
               verbose=verbose-2)

    ##Expr network
    exprNetwork <- WGCNA::adjacency(datExpr=t(genExpr), type="unsigned",
                                    corFnc="cor", corOptions=corInput, power=1)
    message.if(paste("exprNetwork is made, dim:", paste(dim(exprNetwork), collapse="*")),
               verbose=verbose-2)

    ##dnam network
    dnamNetwork <- WGCNA::adjacency(datExpr=eigenloci, type="unsigned",
                                    corFnc="cor", corOptions=corInput, power=1)
    message.if(paste("dnamNetwork is made, dim:", paste(dim(dnamNetwork), collapse="*")),
               verbose=verbose-2)
    message.if("Adjacency matrices were computed.", verbose=verbose-1)

    powerVector <- c()
    outliersNumber <- c()
    result[["corMethod"]] <- corMethod
    result[["mus"]] <- mus
    netsPath <- file.path(outPath,"nets")
    dir.create(netsPath, showWarnings=verbose>2)
    mu2modules <- c()
    for(muValue in mus){
        message.if(paste("mu value : ", muValue), verbose=verbose-1)
        combined <- Pigengene::combine.networks(nets=list(dnamNetwork, exprNetwork),
                                                contributions=c(muValue, 1-muValue),
                                                midfix=muValue, outPath=netsPath,
                                                RsquaredCut=RsquaredCut,
                                                minModuleSize=minModuleSize,
                                                doRemoveTOM=doRemoveTOM,
                                                datExpr=eigenloci, verbose=verbose-3,
                                                doSave=doSaveCombined, doReturNetworks=doReturNetworks)
        powerVector[paste(muValue)] <- combined[["power"]]
        outliersNumber[paste(muValue)] <- length(which(combined$modules==0))
        mu2modules <- rbind(mu2modules, combined$modules)
        rownames(mu2modules)[nrow(mu2modules)] <- muValue
        rm(combined)
        gc()
    }##for
    
    result[["mu2modules"]] <- mu2modules
    ##plot power values per mus
    png(filename=file.path(outPath, paste0("Power-mu-plot.png")))
    plot(x=mus, y=powerVector, xlab="mu", ylab="power")
    dev.off()

    ##plot number of outlier per mus
    png(filename=file.path(outPath, paste0("outliers-mu-plot.png")))
    plot(x=mus, y=outliersNumber, xlab="mu", ylab="Number of outlier genes",
         main=paste("Total number of genes:", nrow(genExpr)))
    dev.off()
    message.if(paste("Two mu plots were saved in:", outPath), verbose=verbose-1)
    
    result <- c(result, powerVector=powerVector, outliersNumber=outliersNumber,
                netsPath=netsPath)
    if(doReturNetworks) ## Add them to result.
        result <- c(result, exprNetwork=exprNetwork, dnamNetwork=dnamNetwork)
    gc()
    return(result)
}
