infer.eigengenes <- function(inetgrator, dnam=NULL, expr, verbose=0){
    ## Takes as input an inetgrator, dnam, and expr data, and computes the eigengenes
    ##corresponding to the inetgrator.
    ## Samples must be on rows.
    
    if(!is.null(dnam)){##if dnam data is provided
        if(any(rownames(dnam)!=rownames(expr)))
            stop("dnam and expr must have the same rownames! ")
    }
    result<-list()
    mu <- inetgrator$mu
    modules <- inetgrator$modules
    eigengeneS <- c() ## Will be a matrix, with samples on rows.
    missingInExpr<-list()
    for(m1Name in names(modules)){
        m1<-modules[[m1Name]]
        message.if(cat("There are",length(names(m1$genes)),"genes in module",m1Name,".\n"), 
                   verbose=verbose)
        method <- colnames(m1$weights)
        
        ## Use expression data?
        if(length(grep(method, pattern="e"))>0){
            DataExpr <- expr[ ,intersect(names(m1$genes), colnames(expr)), drop=FALSE]
            message.if(cat(ncol(DataExpr), "of these genes available in expr.\n"), 
                       verbose=verbose)
            if(ncol(DataExpr)==0)
                warning("No gene of", m1Name, "is available in expr!")
            NAs<-names(m1$genes)[is.na(match(names(m1$genes),colnames(expr)))]
            if(length(NAs)>0){
                DataExpr<-cbind(DataExpr,matrix(0,nrow=nrow(DataExpr),ncol=length(NAs),
                                                dimnames=list(rownames(DataExpr),NAs)))
                message.if(cat("The following genes are missing in expr: \n", NAs, sep="\n"),
                           verbose=verbose-2)
                missingInExpr[[m1Name]] <- NAs
            }
            Data <- DataExpr
        }

        ## Use dnam data?
        if(length(grep(method, pattern="m"))>0){
            dnaM1 <- c() ## The matrix of dnam for the genes in m1 at the gene level.
            message.if("Projecting dnam data...", verbose=verbose)
            er1 <- paste0("DNA methylation data is needed to infer the eigengene for module: ",
                         m1Name, ". Thus, dnam cannot be NULL!")
            if(is.null(dnam))
                stop(er1)
            for(g1 in names(m1$genes)){
                message.if(cat("Projecting ", g1), verbose=verbose-1)
                genes1<-m1$genes[[g1]]
                if(nrow(genes1$loci)==1)
                {
                    toAdd<-dnam[ ,rownames(genes1$loci), drop=FALSE]
                }else{
                    projected <- project.eigen(Data=dnam, pigengene=m1$genes[[g1]]$pigengene, saveFile=NULL)
                    toAdd<-projected$projected
                }
                dnaM1 <- cbind(dnaM1, as.matrix(toAdd))
                colnames(dnaM1)[ncol(dnaM1)] <- g1 
            }
            Data <- dnaM1
        }

        ## Use both dnam and expr:
        if(length(grep(method, pattern="m")>0) & length(grep(method, pattern="e"))>0){            
            Data <- (1-mu)*DataExpr + mu * dnaM1
        }
        
        ## e1 is the inferred eigengene for m1.
        e1 <- project.eigen(Data=Data, pigengene=inetgrator$pigengene[[method]], saveFile=NULL) 
        eigengeneS <- cbind(eigengeneS, as.matrix(e1$projected[,m1Name,drop=FALSE]))
        colnames(eigengeneS)[ncol(eigengeneS)] <- m1Name
    }
    ##End     for(m1Name in names(modules)).
    
    result[["inetgrator"]] <- inetgrator
    result[["eigengeneS"]] <- eigengeneS ## A matrix with samples on rows, the inferred eigengenes on columns.
    result[["missingInExpr"]] <- missingInExpr
    return(result)
}
