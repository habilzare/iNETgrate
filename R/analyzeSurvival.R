analyzeSurvival <- function(survival, event="Dead", Labels, Data=NULL, 
                              excludeSamples=rownames(survival)[as.numeric(survival[, "Time"])<=0],
                              mus, netPath, outPath,
                              favRisk="High", subSet=NULL, doCox=TRUE, 
                              eOrMs="e", time2day=1, until=1,
                              xmax1=max(survival[, "Time"])*(time2day/365), xmin1=0,
                              minRecall4L=0.2, minRecall4H=0.05,
                              verbose=0, doDebug=FALSE){
    ## It performs Cox analysis and AFT (based on e.g., the eigengenes from
    ## the combined network of dnam and expr data)
    ## Inputs:
    ## survival: A data frame with patients on rows.
    ##^Columns include: "Time" (in days if time2day=1), event (e.g., "Dead"), and some risk criteria.
    ##Abnormality and Risk columns (columns that start with Risk e.g. Risk1, Risk2, ...) are
    ##optional.
    ## time2day: Set it to 1 if the "Time" column in survival is in days. Otherwise, this is the
    ## number of days in a unit of "Time". E.g., set it to 30 if Time is in months.  
    ## mus: A vector of mu values.
    ## feature: A matrix contains eigengenes.
    ## netPath: A string to showed the path of eigengenes. If NULL, the Data matrix is used.
    ## Data: The feature matrix, e.g., patients on rows and eigengenes on columns. Exactly one of
    ## netPath or Data must be not NULL.
    ## outPath: A string to showed the path for saved object.
    ## favRisk: A risk category which is predicted by the network classifier and
    ##is going to be analized more precisely. The value should be in c("High", "Int", "Low").
    ## subSet: A risk category which is predicted by the other classifiers and
    ##is going to be analyzed more precisely. The value should be in c("High", "Int", "Low").
    ## doCox: A boolean flag which is by defualt true and makes the function to perform cox
    ## analysis.
    ## eOrMs: A character vector that determines what data were used to compute eigengenes,
    ##^e.g., c("m", "e", "em"). It will be ignored if netPath is NULL.
    ## Labels: A vector with only "High", "Int", and "Low" values. Names must be similar to the row
    ## names of survival.
    ## Output:
    ## bestPvalues:
    ## bestModules:
    ## kmPVals:
    ## until: See accelFailAnalysis().
    starTime <- Sys.time()
    if(is.null(netPath) && is.null(Data))
        stop("Both netPath and Data can not be NULL !")
    if(!is.null(netPath) && !is.null(Data))
        stop("One of the netPath or Data should be NULL !")
    if(!is.null(favRisk) && !(favRisk %in% c("High", "Int", "Low")))
        stop("The value of favRisk should be in High, Int or Low")
    if(!is.null(subSet) && !(subSet %in% c("High", "Int", "Low")))
        stop("The value of subSet should be in High, Int or Low")
    class(survival[, "Time"]) <- "numeric"
    class(survival[, event]) <- "numeric"

    ## Accelerated failure time and Cox analyses on all data
    message.if("Accelerated failure time and Cox analysis ...", verbose=verbose)

    ## defining input time(years or days?!) as an input for coxAnalysis
    inputTime <- survival[, "Time"] * time2day
    names(inputTime) <- rownames(survival)

    ## omit those patients with inputTime<=0
    keepSamples <- names(inputTime[inputTime>0])
    message.if(verbose=verbose-1,
               paste("The survival time is 0 for the following", sum(inputTime==0),
                     "patients. Excluded."))
    ## Exclude some samples?
    keepSamples <- keepSamples[!keepSamples %in% excludeSamples]
    message.if(verbose=verbose-1, names(inputTime)[inputTime==0])
    inputTime <- inputTime[keepSamples]
    inputTime <- inputTime/365
    inputEvent <- survival[keepSamples, event]
    names(inputEvent) <- keepSamples
    message.if(paste("inputEvent: \n"), verbose=verbose)
    if(verbose)
        print(table(inputEvent))

    afts <- list()
    bestPvalues <- data.frame()
    returnbestModules <- list()
    returnkmPVals <- list()

    aftPaths <- c() ## A matrix with mu on rows and risk on columns.
    pastedeOrMs <- paste(eOrMs, collapse="-")
    for(mu in mus){
        muChar <- as.character(mu)
        message.if(paste0("mu value : ", mu), verbose=verbose)
        if(!is.null(netPath)){ ## Load the eigengenes and cbind them if needed.
            muPath <- file.path(netPath, paste0("mu", mu))
            coxData <- c()##cbind all the above matrices. "e", "m", "em".
            for(eOrM in eOrMs){
                fileIn <- file.path(muPath,  paste0("Pigen_", eOrM), "eigengenes.RData")
                eigengenes <- get(load(file=fileIn, verbose=TRUE)) ##eigengenes
                ##^each file has a matrix of eigengenes for dnam, expr, dnam&expr
                coxData <- cbind(coxData, eigengenes)
                rm(eigengenes)
            }
        } else { ## Do not load eigengenes,
            coxData <- Data
        }
        message.if(verbose=verbose-2, "Dim of original coxData:")
        message.if(verbose=verbose-2, paste(dim(coxData), collapse="*"))
        ## Rearrange input matrix as an input for coxAnalysis
        coxData <- coxData[rownames(coxData) %in% keepSamples, , drop=FALSE]
        message.if(verbose=verbose-2, paste("Dim of original coxData after keeping only ",
                                            length(keepSamples),
                                            "samples with >0 time data."))
        message.if(verbose=verbose-2, paste(dim(coxData), collapse="*"))
        ## QC:
        if(nrow(coxData)==0){
            stop("Check row names of survival and the (Data) eigengenes matrix, do they overlap?!")
        }

        toKeepCol <- c()## we don't like all NA coloumns
        for(i1 in seq_len(ncol(coxData))){
            if(sum(is.na(coxData[, i1]))!=nrow(coxData))
                toKeepCol <- c(toKeepCol, i1)
        }
        message.if(verbose=verbose-2, paste(ncol(coxData)-length(toKeepCol),
                                            "eigengenes were all NA, and removed"))        
        coxData<-coxData[, toKeepCol, drop=FALSE]

        ##cox
        if(doCox){
            bestModules <- coxAnalysis(dataCategory="allSamples",
                                       inputTime=inputTime[rownames(coxData)],
                                       inputEvent=inputEvent[rownames(coxData)],
                                       inputMatrix=coxData,
                                       outPath=outPath, verbose=verbose-1)
            returnbestModules[[paste0("mu", mu)]] <- bestModules
        } else {
            bestModules <- colnames(coxData)
        }
        inputTime <- inputTime[rownames(coxData)]
        inputEvent <- inputEvent[rownames(coxData)]
        ## ###################################################
        ## ####### Accelerated Failure Time analysis #########
        ##  ##################################################
        selectedModules <- bestModules[seq_len(min(3, length(bestModules)))]
        sampleData <- data.frame(coxData[, selectedModules, drop=FALSE], "time"=inputTime, "event"=inputEvent)
        names(sampleData) <- make.names(names(sampleData))
       
        yAllSamples <- Surv(time=inputTime, event=inputEvent, type="right")
        right1 <- paste(make.names(colnames(coxData[, selectedModules, drop=FALSE])), collapse="+")
        formulaI <- as.formula(paste("yAllSamples ~", right1))

        if(verbose>=3)
            print(formulaI)

        aft <- list()        
        message.if(verbose=verbose-2, paste0("Accelerated Failure Time analysis..."))
        selectedModules <- make.names(selectedModules, unique=TRUE)
        inputEvent <- setNames( as.logical(inputEvent), names(inputEvent))
        ## Change "ME0-e" to "ME0.e" not to cause issues in formulae 
        colnames(coxData) <- make.names(colnames(coxData), unique=TRUE)

        ## KM plot of all samples:
        message.if("KM plot of all the samples:", verbose=verbose-1)
        yTrain <- Surv(time=inputTime, event=inputEvent, type="right")
        names(yTrain) <- names(inputTime)
        condA <- rep(x="all",times=length(inputTime))
        png(file.path(outPath, "allSamplesKM.png"))        
        survfitOutputAll <- survfit(formula=yTrain ~ condA, data=as.data.frame(condA))
        plot(survfitOutputAll, xlim=c(xmin1, min(xmax1, max(inputTime, na.rm=TRUE))),
             col='blue', main=NULL, xlab="Time (years)", ylab="Survival Probability")
        dev.off()

        message.if(verbose=verbose-2, paste("AFT based on survival time and", event))
        riskLevel <- c("High"="High","Low"="Low","Int"="Int")
        aftPathsL <- file.path(outPath, paste0("mu", mu),
                               paste0("survival_", pastedeOrMs))

        ## QC:
        if(!all(unique(Labels) %in% c("High", "Int", "Low")))
           stop("Labels values must be High, Int, or Low!")
        if(any(! rownames(coxData) %in% names(Labels)))
           stop("Row names of Data must be a subset of names of Labels. Some samples have no Labels!")
        survivalCox <- cbind(survival[rownames(coxData), , drop=FALSE], Risk=Labels[rownames(coxData)])
        if(doDebug){
            message.if(verbose=verbose-1, "Run AFT?")
            
        }
        aft <- accelFailAnalysis(survival=survivalCox, eventCol=event,
                                 time2day=time2day,
                                 Data=coxData[ ,selectedModules, drop=FALSE], xmax1=xmax1, 
                                 minRecall4L=minRecall4L, minRecall4H=minRecall4H, until=until,
                                 resultPath=aftPathsL,
                                 riskCol="Risk",
                                 abnormalityCol=NULL,
                                 doAddConfTable=TRUE, xmin1=xmin1,
                                 risk2col=c( "Low"="green", "Int"="blue", "High"="red"),
                                 favRisk=favRisk, riskLevel=riskLevel,
                                 subSet=subSet,
                                 ylabKm="Survival probability", verbose=verbose-2)
        ## ## ######### Accelerated Failure Time Done ############
        ##best p-values for the mu value:
        kmPVals <- aft$kmPVals
        returnkmPVals <- as.matrix(cbind(returnkmPVals, kmPVals))
        colnames(returnkmPVals)[ncol(returnkmPVals)] <- muChar
        bestP <- min(kmPVals)
        bestFeaturesP <- names(sort(kmPVals))[1]
        toAdd <- data.frame("pvalue"=bestP, "mu"=mu, "features"=bestFeaturesP)
        bestPvalues <- rbind(bestPvalues, toAdd)
        aftPaths <- c(aftPaths, aftPathsL)
        names(aftPaths)[length(aftPaths)] <- muChar
        afts[[muChar]] <- aft
    }
    ##^ for(mu in mus)
    
    ##plot
    bPFileName <- paste0("bestPvalues_", pastedeOrMs)
    png(filename=file.path(outPath, paste0(bPFileName, ".png")),
        height=2*480, width=2*480, res=2.5*72)
    bestPvaluePlot <- plot(x=mus, y=log10(bestPvalues[, "pvalue"]),
                           xlab=expression(mu), ylab="p-value (log10)", cex.lab=1.4)
    dev.off()

    ## Conclusion:
    write.csv(bestPvalues, file=file.path(outPath, paste0(bPFileName, ".csv")))
    if(verbose>1)
        print(bestPvalues)
    
    ## Timing:
    timeTaken <- Sys.time() - starTime
    message.if(me=paste("Survival analysis took:", timeTaken, attr(timeTaken,"units"), "\n"), 
               verbose=verbose-3) 
    ## Output:
    return(list(mus=mus,
                bestPvalues=bestPvalues,
                bestModules=returnbestModules,
                kmPVals=returnkmPVals,
                bestPvaluePlot=bestPvaluePlot,
                afts=afts,
                aftPaths=aftPaths,
                timeTaken=timeTaken))
    }
