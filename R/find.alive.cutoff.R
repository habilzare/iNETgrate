find.alive.cutoff <- function (hope, time1, until, minRecall=0.2, risk="Low", Labels=NULL, 
                               doDebug=FALSE, resPath, verbose=0, doPlot=FALSE,
                               lowLabel="Low", highLabel="High"){
    ## time1: A numeric vector.
    ## hope: A numeric vector index the same as time1.
    ## ratio: A numeric in the range [0, 1]. Depreciated after minRecall was added.
    ## until: A numeric (10 years for breast cancer, 2 years for AML). It will 
    ## be ignored if Labels is not null.
    ## lowLabel: The label of low risk-cases, e.g., "Favorable".
    ## highLabel: The label of  high-risk cases, e.g., "Poor".
    ## Labels : A vector named the same as time1 and hope, which has the training risk groups:
    ## "Poor", "Intermediate/Normal", "Favorable"

    ## risk: Character. Set to "Low" or "High".
    ## Assumption: The number of alive cases with higher hope values should be generally more.
    ## Finds the best (smallest) cutoff on hope where the ratio of alive cases, for which time1 > until, 
    ##^is at least 'ratio'.
    ## resPath: Path to the folder to save results (plots)  
    result <- list()
    message.if("Finding alive cutoff...", verbose=verbose)
    message.if(paste("risk:", risk), verbose=verbose-1)
    result[["call"]] <- match.call()
    if(any(names(hope) != names(time1)))
        stop("Names of hope and time1 must be identical in the same order!")
    o1 <- order(hope)
    hope <- hope[o1]; time1 <- time1[o1] ##; dead <- dead[o1]
    highPrecision <- c()
    lowPrecision <- c()
    highRecall <- c()
    lowRecall <- c()
    ind0 <- 0
    maxHighPrecision <- 0
    maxLowPrecision <- 0
    num <- NULL

    ##QC
    if(!is.null(until)&!is.null(Labels))
        warning("Because Labels is not NULL, until is ignored!")
    if(is.null(until)&is.null(Labels))
        stop("Both Labels and until cannot be NULL!")
    if(is.null(resPath))
        stop("resPath is null!! Choose a folder to save results")
    for(ind in seq_along(hope)){
        cutoff <- hope[ind]
        message.if(paste("cutoff:", cutoff), verbose=verbose-3)
        if(is.null(Labels)){
          tpHigh <- sum(time1[hope<=cutoff]<until)
          tpLow <- sum(time1[hope>=cutoff]>until)
          pHigh <- sum(time1<=until)
          pLow <- sum(time1>until)
          if(pLow==0)
              stop("There is no value in time1 that is greater than until!")
          if(pHigh==0)
              stop("There is no value in time1 that is smaller than until!")
        }else{ ## Labels is not NULL,
          Labels<-Labels[names(hope)]
          tpHigh <- sum(hope<=cutoff & Labels==highLabel)
          tpLow <- sum(hope>=cutoff & Labels==lowLabel)
          pHigh <- sum(Labels==highLabel)
          pLow <- sum(Labels==lowLabel)
          if(pLow==0)
              stop("Labels has no value equal to lowLabel!")
          if(pHigh==0)
              stop("Labels has no value equal to highLabel!")
        }
        highPrecision[ind] <- tpHigh/sum(hope<=cutoff)
        lowPrecision[ind] <- tpLow/sum(hope>=cutoff)
        highRecall[ind] <- tpHigh/pHigh
        lowRecall[ind] <- tpLow/pLow
        
        if(risk=="High"){ ## High-risk
            if(highRecall[ind] >= minRecall & highPrecision[ind] >= maxHighPrecision){
                ind0 <- ind
                maxHighPrecision <- highPrecision[ind]
                num <- ind0
            }
        }else{ ## Low-risk
            if(lowRecall[ind] >= minRecall & lowPrecision[ind] >= maxLowPrecision){
                ind0 <- ind
                maxLowPrecision <- lowPrecision[ind]
                num <- length(hope)-ind0+1
            }
        }
        if(doDebug){
            print(paste("cutoff:", cutoff))
            print(paste("tpLow:", tpLow))
            print(paste("lowPrecision[ind]:", lowPrecision[ind]))
            plot(x=hope[seq_len(ind)], y=time1[seq_len(ind)], xlim=range(hope), ylim=range(time1))
            title(paste("maxLowPrecision=", maxLowPrecision, ", maxHighPrecision=", maxHighPrecision))
            browser()
        }
    }
    cutoff <- hope[ind0]
    if(doPlot){
        if(risk=="Low"){
            plotFile <- file.path(resPath, "precision_recall_lowRisk.png")
            png(filename=plotFile)
            plot(lowRecall, lowPrecision)
            points(lowRecall[ind0], lowPrecision[ind0], col='red')
            dev.off()
        }
        if(risk=="High"){
            plotFile <- file.path(resPath, "precision_recall_highRisk.png")
            png(filename=plotFile)
            plot(highRecall, highPrecision)
            points(highRecall[ind0], highPrecision[ind0], col='red')
            dev.off()
        }
        message.if(paste("Plotted: ", plotFile), verbose=verbose-1)
        result[["plotFile"]] <- plotFile
    }

    ##
    message.if(paste("cutoff:", cutoff), verbose=verbose-1)
    result[["num"]] <- num
    result[["cutoff"]] <- cutoff
    result[["highPrecision"]] <- highPrecision
    result[["lowPrecision"]] <- lowPrecision
    result[["highRecall"]] <- highRecall
    result[["lowRecall"]] <- lowRecall
    return(result)
}
