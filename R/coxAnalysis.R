coxAnalysis <- function(dataCategory="training", inputTime,
                         inputEvent, inputMatrix, outPath, weight=NULL,
                         verbose=0){
    message.if("Cox analysis...", verbose=verbose)
    if(inherits(inputMatrix, "data.frame")){
        inputMatrix <- as.matrix(inputMatrix)
    }
    if(!inherits(inputMatrix, "matrix")){
        stop("inputMatrix must be a matrix or a data.frame!")
    }

    ## Output Files
    message.if(paste("coxAnalysis results will be saved at:", 
                     outPath, "\n"), verbose=verbose-1)
    glmAllFile <- file.path(outPath, paste0(dataCategory, "_coefGlmnet.txt"))
    
    y1Out <- file.path(outPath, paste0(dataCategory, "_y.RData"))
    
    pngGlmOut <- file.path(outPath, paste0(dataCategory, "_regularization.png"))

    dir.create(path=outPath, recursive=TRUE, showWarnings=FALSE)
    
    ## Calculate Surv object. This is dependent on the category of data. Also
    ## create a temporary yVar variable
    if(dataCategory=="training"){
        yTrain <- Surv(time=inputTime, event=inputEvent, type="right")
        save(yTrain, file=y1Out)
        yVar <- yTrain
    } else if(dataCategory=="validation"){
        yValid <- Surv(time=inputTime, event=inputEvent, type="right")
        save(yValid, file=y1Out)
        yVar <- yValid
    } else if (dataCategory=="allSamples") {
        yAllSamples <- Surv(time=inputTime, event=inputEvent, type="right")
        save(yAllSamples, file=y1Out)
        yVar <- yAllSamples
    }

    ## Weight of 1 to all if NULL
    if(is.null(weight)){
        weight <- c(rep(1,nrow(yVar)))
    }
    ## Penalized Cox
    fitGlmnet <- glmnet::glmnet(x=inputMatrix, y=yVar, family="cox", alpha=1,
                                weights=weight);

    ## coef
    capture.output(coef(fitGlmnet)[,seq_len(10)], cat("\n\n"),
                   file=glmAllFile, append=TRUE)
    
    png(pngGlmOut)
    plot(fitGlmnet)
    dev.off()

    ## Necessary variable for upcoming logic
    coGlmMatrix <- coef(fitGlmnet)

    ## Find the best modules, settings ::
    ## Determine best modules
    ind <- 1
    
    ## Index for bestModules
    bestModIndex <- 1
    
    ## Vector to hold the best modules
    bestModules <- character(nrow(coGlmMatrix))
    orderedCoeff <- c()

    while((ind < ncol(coef(fitGlmnet))) & !is.vector(coGlmMatrix)){
        sortedCoef <- sort(abs(coGlmMatrix[,ind]), decreasing=TRUE)
        sortedCoef <- sortedCoef[!(sortedCoef==0)]
        while(length(sortedCoef) !=0){
            bestModules[bestModIndex] <- names(sortedCoef[1])
            orderedCoeff[names(sortedCoef[1])] <- sortedCoef[1]
            coGlmMatrix <- coGlmMatrix[!rownames(coGlmMatrix) %in%
                                       names(sortedCoef[1]),];
            
            sortedCoef <- sortedCoef[!names(sortedCoef) %in%
                                     names(sortedCoef[1])];
            bestModIndex <- bestModIndex + 1
            if(is.vector(coGlmMatrix)){
                break
            }
        }
        ind <- ind + 1
    }
    ## Get last module
    bestModules[bestModules==""] <- rownames(coef(fitGlmnet))[!rownames(coef(fitGlmnet)) %in%
                                                              bestModules]
    save.if(orderedCoeff, file=file.path(outPath, paste0(dataCategory,"_coeff.RData")),
            verbose=verbose-1)    
    
    return(bestModules)
}
