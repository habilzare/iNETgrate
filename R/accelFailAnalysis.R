accelFailAnalysis <- function(Data, survival, time2day, eventCol="Dead",
                              riskCol="Risk1", weight=NULL, minRecall4L=0.2,
                              minRecall4H=0.1, until=10, xmin1=0,
                              xmax1=max(survival[, "Time"])*(time2day/365),
                              predType="lp", doTitle=TRUE, resultPath, 
                              risk2col=c("Low"='green', "Int"='blue', "High"='red'),
                              doAddConfTable=FALSE, 
                              favRisk="High", riskLevel, subSet=NULL, pvalDigits=0, 
                              ylabKm="Survival probability", verbose=0){

    ## Performs accelerated failure time model and prediction on each combination of
    ##^selected modules.
    ## inputParam is the set with the best modules
    ## trainParam is a list with the names of the training dataset
    ## survival includes the time as well as the
    ## minRecall4L is our cutoff recal for low risk predictions
    ## minRecall4H is our minimum recal for high risk predictions
    ## until is the minimum survival time ***in years*** for an alive case to be considered as low risk.
    ## set it as NULL if you want to create km-plots based on riskCol
    ## xmax1 is the maximum time for plotting the K-M curves in years.
    ## Data are the (projected) eigengenes for each sample, including the projected values

    ##riskLevel: an ordered named vector which values are risk group names(subset of unique(riskCol)) and
    ##names are Low, Int and High.
    ##This is neccesary for  computing confusion matrix
    ##riskCol: the col in the survival that has the risk groups, the cytogenetic col for the AML survival study
    ##favRisk : the risk group that you want to draw kmplot for their real risk level.
    ##It should be Low, High or Int.
    ## subSet: A character the determines which risk category from the survival data will be further
    ##stratified based on the predicted risk. E.g., "Int" leads to plotting the predicted risk
    ##groups for the subset of data that are intermediate risk based on survival data.  
    ##inputTime: inputTime in years.
    ##survival: must have "Time" col.
    ##yTrain: Used to be an input, now deprecated.
    ## pvalDigits: Number of digits of the p-value shown on the KM plots.
    ## Set it to NA not to add pvalues on the plots.
    
    result <- list()
    message.if("Accelerated failure time Analysis ...", verbose=verbose)

    ## QC:
    if(!"Time" %in% colnames(survival))
        stop("survival should have a 'Time' column!")
    if(!eventCol %in%colnames(survival))
        stop("survival should have a ", eventCol, " column!")
    if(!riskCol %in% colnames(survival))
        stop("survival should have a ", riskCol, " column!")
    ## Checking rows of survival and Data!
    if(nrow(survival)!=nrow(Data)){
        inBoth <- intersect(rownames(survival), rownames(Data))
        if(length(inBoth)==0)
            stop("Row names in Data and survival should intersect!")
        message.if(paste("Analyzing only ", length(inBoth), "samples for which \n",
                       "both survival and eigengene data are available.\n"),
                   verbose=verbose-1)
        survival <- survival[inBoth, ]
        Data <- Data[inBoth, ]
    }

    inputTime <- as.numeric(survival[, "Time"]) * time2day
    ## inputEvent <- as.numeric(survival[, eventCol])
    names(inputTime) <- rownames(survival)

    ## omit those patients with inputTime<=0
    keepSamples <- names(inputTime[inputTime>0])
    message.if(verbose=verbose-1,
               paste("The survival time is 0 for the following", sum(inputTime==0),
                     "patients. Excluded."))
    ## Exclude some samples?
    ## keepSamples <- keepSamples[!keepSamples %in% excludeSamples]
    message.if(verbose=verbose-1, names(inputTime)[inputTime==0])
    inputTime <- inputTime[keepSamples]
    inputTime <- inputTime/365
    inputEvent <- survival[keepSamples, eventCol]
    names(inputEvent) <- keepSamples
    message.if("inputEvent: ", verbose=verbose)
    ##eventMat <- t(as.matrix(table(inputEvent)))
    ##message.if(me=paste(colnames(eventMat), eventMat), verbose=verbose)
    message.if(print(table(inputEvent)), verbose=verbose)
    survival <- survival[keepSamples, ]

    ##QC:
    inputEvent <- inputEvent[names(inputTime)]
    if(length(inputTime)!= length(inputEvent))
        stop("Do names match in inputTime and inputEvent?!")

    ## Change "ME0-e" to "ME0.e" not to cause issues in formulae 
    colnames(Data) <- make.names(colnames(Data), unique=TRUE)
    ## riskCol <- make.names(riskCol, unique=TRUE)
    ## colnames(survival) <- make.names(colnames(survival), unique=TRUE)

    message.if(me=paste("Risk levels", riskLevel), verbose=verbose-2)
    if(!("Int" %in% names(riskLevel)))
        riskLevel <- c(riskLevel["Low"], "Int"="Int", riskLevel["High"])
    ##^e.g., riskLevel["High"]="Poor"
    levelRisk <- names(riskLevel)
    names(levelRisk) <- riskLevel
    ##^e.g., levelRisk["Poor"]="High"
    if(favRisk!="High" & favRisk!="Low" & favRisk!="Int")
        stop("favRisk should be Low, High or Int")
    yTrain <- Surv(time=inputTime, event=inputEvent, type="right")
    names(yTrain) <- names(inputTime)
    dir.create(path=resultPath, recursive=TRUE, showWarnings=FALSE)
    message.if(paste("accelFailAnalysis results will be saved at:", 
                     resultPath), verbose=verbose-1)


    ## Int risk patients based on cytogenetic 
    subsetPatients <- NULL
    if(!is.null(subSet)){
        subsetPatients <- rownames(survival)[which(survival[, riskCol]==subSet)]
    }

    ## If weight is not specified
    if(is.null(weight)){
        weight <- inputTime ## To match the names
        weight[] <- 1
    }

    ## Create Power Set
    powerSet <- list()
    featureNames <- colnames(Data)
    for(counter in seq_along(featureNames)) {
        powerSet <- c(powerSet, combn(featureNames, m=counter, simplify=FALSE))
    }

    ## Useful survival data: 
    ## trainSurvival <- survival[rownames(survival) %in% rownames(Data), ]
    trainSurvival <- cbind(survival, "Years"=survival[, "Time"]*time2day / 365)
    trainSurvival[, eventCol] <- (trainSurvival[, eventCol] == 1) * 1
    ## Hanie : Risk is also defiend for trainsurvival

    ## Concatenate
    ## Note: Subsetting orders the data to eigengenes, which is necessary
    trainD <- as.data.frame(cbind(Data[rownames(trainSurvival), , drop=FALSE], trainSurvival))

    ## Fit accelerated failure time model to all sets and plot
    message.if(me="Fitting accelerated failure models....", verbose=verbose)

    ## save the p-values in this object
    kmPVals <- c()
    subsets <- list()
    for(currentModule in powerSet){
        ## Rename the element of power set
        subPowerSet <- paste(currentModule, collapse="_")
        message.if(paste("Working on:", subPowerSet), verbose=verbose-1)
        subsets[[subPowerSet]] <- list()
        kmPVals[subPowerSet] <- 1
        
        ## Output .txt files
        subsetPath <- file.path(resultPath, subPowerSet)
        dir.create(path=subsetPath, recursive=TRUE, showWarnings=FALSE)
        outFile <- file.path(subsetPath, paste("cox_aft.txt", sep=""))
        message.if(paste("accelFailAnalysis is reporting in: ", outFile), verbose=verbose-1) 

        ## Surv
        formulaI <- as.formula(paste("yTrain ~", paste(currentModule, collapse="+")))
        if(verbose >=3)
            print(formulaI)
        fitI <- survreg(formula=formulaI, 
                        data=as.data.frame(Data[rownames(survival), , drop=FALSE]), 
                        dist='weibull', scale=1, weights=weight)

        capture.output(summary(fitI), file=outFile)
        message.if(paste("A summary of the survival model was written in", outFile), 
                   verbose=verbose-1)
        currentModulePasted <- paste(currentModule, collapse="_")
        ## save model for future analysis
        subsets[[subPowerSet]][["fitI"]] <- fitI

        
        ## Make predictions and find good cutoffs
        predictI <- predict(fitI, type=predType, newdata=as.data.frame(trainD))
        if(any(is.na(predictI)))
            stop("NA predicted for a case!")
        doRisk <- FALSE
        if(is.null(until)){
            riskLabels <- as.character(survival[, riskCol])
            names(riskLabels) <- rownames(survival)
            doRisk <- TRUE
        }

        message.if("Finding cutoffs...", verbose=verbose-1)
        customized.fac <- function(custMinRecall, custRisk, custLabels=NULL){
            time1 <- setNames(trainD[, "Years"], nm=rownames(trainD))
            f1 <- findAliveCutoff(hope=predictI, time1=time1, 
                                    until=until, verbose=verbose-3, resPath=subsetPath, 
                                    Labels=custLabels, 
                                    minRecall=custMinRecall, risk=custRisk)
            return(f1)
        }
        if(doRisk){
            found <- customized.fac(custLabels=riskLabels, custMinRecall=minRecall4L, custRisk="low")
            foundH <- customized.fac(custLabels=riskLabels, custMinRecall=minRecall4H, custRisk="High")
        }else{
            found <- customized.fac(custMinRecall=minRecall4L, custRisk="low")
            foundH <- customized.fac(custMinRecall=minRecall4H, custRisk="High")
        }

        ## Keep our cutoffs for future analysis
        subsets[[subPowerSet]][["found"]] <- found
        subsets[[subPowerSet]][["foundH"]] <- foundH
        
        
        ## Condition for training:
        ## cond is a vector of "High", "Low", and "Int" values named by the cases names.
        cond <- rep(x=NA, times=length(predictI))
        cond[predictI>=found$cutoff] <- "Low"
        cond[predictI<=foundH$cutoff] <- "High"
        cond[is.na(cond)] <- "Int"
        if(found$cutoff < foundH$cutoff){
            cond[predictI>=found$cutoff & predictI<=foundH$cutoff] <- "Int"
            warning("The lower cutoff is unusually smaller than the higher cutoff.", 
                    " The cases in the overlap were considered Int risk!")
        }
        names(cond) <- names(predictI)
        if(any(rownames(trainD) != names(cond))) {
            stop("Rows of trainD are assumed to be in the same order as names of cond!")
        }
        trainDset <- cbind(trainD, cond=cond[rownames(trainD)])
        subsets[[subPowerSet]][["fitI"]] <- fitI

        outCondMat <- t(as.matrix(table(cond)))
        message.if(me=paste(colnames(outCondMat), outCondMat),
                   verbose=verbose-2)

        ## Did we get more than 1 condition?
        nonIntPatients <- names(which(cond!="Int"))
        if(length(unique(cond[nonIntPatients]))<2){
            warning(subPowerSet, 
                    ": The cases were categorized into only one condition!? Skipped!")
            next
        } else {
            message.if("Log-Rank test to calculate P-value for training dataset...", 
                       verbose=verbose-1)
            ## We don't consider Medium (Int) risk in computing the p-value.
            survRes <- Surv(time=trainDset[nonIntPatients, "Years"], 
                            event=trainDset[nonIntPatients, eventCol], type='right')
            survdiffRes <- survdiff(formula=survRes ~ cond[nonIntPatients])
            subsets[[subPowerSet]][["survdiffRes"]] <- survdiffRes
            capture.output(survdiffRes, file=outFile, append=TRUE)
        }
        
        ## Create scatter plots showing predicted time vs. actual time, colored by the event.
        scatterFile <- file.path(subsetPath, "comparison.png")
        png(scatterFile)
        colorChart <- trainDset[, eventCol]
        colorChart[which(colorChart==1)] <- 2 ##red for died
        colorChart[which(colorChart==0)] <- 3 ##green for alive
        plot(predictI, trainDset[, "Years"], col=colorChart, 
             xlab="Predicted Time (Years)", ylab="Time To Last Followup (Years)", cex.lab=1.4)
        abline(v=found$cutoff, col='darkgreen')
        abline(v= foundH$cutoff, col='orange')
        titleL <- paste("# of predicted low-risk cases (recall >", minRecall4L, ")=", found$num)

        titleH <- paste("# of predicted high-risk cases (recall >", minRecall4H, ")=", foundH$num)
        if (doTitle) {
            title(paste(titleL, titleH, sep="\n"))
        }
        legend("topright", col=c(2, 3), c("Died", "Alive"), pch=1)
        dev.off()
        message.if(paste("Predicted time vs. actual time were plotted in:", scatterFile), 
                   verbose=verbose-1)

        ## KM Plot based on predictions and also the input groups:
        png(file.path(subsetPath, "comprehensive.png"), width =2800, height=1800, pointsize=20)
        if(doAddConfTable){ ## Do we have the risk information?
            par(mfrow=c(3, 3))
        }
        ## KM plot, predictions for all patients
        plottedAll <- plotKM(inputTime=inputTime, inputEvent=inputEvent, 
                             cond=cond[names(inputEvent)],  
                             titleText=paste("Module", subPowerSet, ": all patients based on iNETgrate"), 
                             cond2Color=risk2col, weight=weight, 
                             plotFile=file.path(subsetPath, "predicted_all.png"), ylab1="", 
                             pvalDigits=pvalDigits, xmax1=xmax1, xmin1=xmin1, verbose=verbose-2)
        kmPVals[subPowerSet] <- plottedAll$pVal
        titleL <- paste("# of predicted low-risk cases (recall >", minRecall4L, ")=", found$num)
        titleH <- paste("# of predicted high-risk cases (recall >", minRecall4H, ")=", foundH$num)
        if (FALSE) { ## Add the above numbers?
            title(paste(titleL, titleH, sep="\n"))
        }
        if(doAddConfTable){
            message.if("Adding the confusion matrix...", verbose=verbose-1)
            ##To add kmplot based on only given risk (e.g., cytogenetics)
            toKeep <- rownames(survival[which(survival[, riskCol]%in%riskLevel), ])
            survival1 <- survival[toKeep, ]
            patientNames <- intersect(names(cond), rownames(survival1))
            survival1 <- survival1[patientNames, ]
            ##
            ## KM plots based on the input risk
            risk1 <- levelRisk[as.character(survival1[patientNames, riskCol])]
            names(risk1) <- patientNames
            plottedIn <- plotKM(inputTime=inputTime[patientNames], inputEvent=inputEvent[patientNames], 
                                cond=risk1, 
                                titleText=paste("All patients based on Labels"), 
                                cond2Color=risk2col, weight=weight[patientNames], 
                                plotFile=file.path(subsetPath, paste0("Labels_all.png")), ylab1=ylabKm, 
                                pvalDigits=pvalDigits, xmax1=xmax1, xmin1=xmin1, verbose=verbose-2)

            ## Confusion matrix:
            rSum <- c()
            cSum <- c()
            cond1 <- cond[patientNames]
            truth <- factor(survival1[, riskCol], labels=names(riskLevel[c("Low", "Int", "High")]), 
                            levels=riskLevel[c("Low", "Int", "High")])
            pred <- factor(cond1, labels=c("Low", "Int", "High"), levels=c("Low", "Int", "High"))
            xtab <- table(pred, truth)

            cfMatrix <- caret::confusionMatrix(xtab)
            cft <- as.matrix(cfMatrix)
            for(i1 in seq_len(nrow(cft))){
                rSum[i1] <- sum(cft[i1, ])
                cSum[i1] <- sum(cft[, i1])
            }
            sumAll <- sum(rSum)
            cSum <- c(cSum, sumAll)
            cft <- cbind(cft, "rowSum"=rSum)
            cft <- rbind(cft, "colSum"=cSum)
            textplot(cft, show.rownames=TRUE, show.colnames=TRUE)
            title("confusion matrix")
            ##Keep cft
            subsets[[subPowerSet]][["confusion"]] <- cft

            ##favRisk : creates a km-plot for patients predicted as favRisk
            favRiskPredict <- names(which(cond==favRisk))
            toKeep <- rownames(survival[which(survival[, riskCol]%in%riskLevel), ])
            survival2 <- survival[toKeep, ]
            patientNames <- intersect(favRiskPredict, rownames(survival2))
            message.if(paste("KM plots for the favorite predicted risk group:", favRisk), verbose=verbose-2) 
            risk2 <- levelRisk[as.character(survival2[patientNames, riskCol])]
            names(risk2) <- patientNames
            pFile2 <- file.path(subsetPath, paste0("Labels_of_predicted_", favRisk, ".png"))
            titleFav <- paste0(length(patientNames), " patients predicted by iNETgrate as ", favRisk, 
                               "-risk classified based on Labels")
            plottedPred <- plotKM(inputTime=inputTime[patientNames], 
                                  inputEvent=inputEvent[patientNames], 
                                  cond=risk2, titleText=titleFav, verbose=verbose-2, 
                                  cond2Color=risk2col, weight=weight[patientNames], 
                                  plotFile=pFile2, ylab1="", pvalDigits=pvalDigits, xmax1=xmax1, xmin1=xmin1)
            ##
            realCond <- survival2[patientNames, ]
            realLow <- rownames(realCond[which(realCond[, riskCol]==riskLevel["Low"]), ])
            if("Int" %in% names(riskLevel))
                realInt <- rownames(realCond[which(realCond[, riskCol]==riskLevel["Int"]), ])
            realHigh <- rownames(realCond[which(realCond[, riskCol]==riskLevel["High"]), ])

            ## KM Plot for a subset of data (e.g., KM plots, for the cytogenetically med cases.)
            if(!is.null(subsetPatients)){
                titleText1 <- paste0(length(subsetPatients), " patients with Labels ", subSet, 
                                     "-risk classified based on iNETgrate")
                plotFile1 <- file.path(subsetPath, paste0("predicted_of_Labels_", 
                                                          subSet, "_subset.png"))
                if(any(!subsetPatients %in% names(inputTime)))
                    stop("subsetPatients must be a subset of names of inputTime!")
                plottedSub <- plotKM(inputTime=inputTime[subsetPatients], 
                                     inputEvent=inputEvent[subsetPatients], 
                                     cond=cond[subsetPatients], 
                                     doIntLow4pval=TRUE, 
                                     titleText=titleText1, 
                                     cond2Color=risk2col, 
                                     weight=weight[subsetPatients], 
                                     plotFile=plotFile1, 
                                     ylab1=ylabKm, pvalDigits=pvalDigits, xmax1=xmax1, xmin1=xmin1, 
                                     verbose=verbose-2)
                survfitOutput <- plottedSub$survfitOutput
                subsets[[subPowerSet]][["survfitOutput"]] <- survfitOutput
            }

            ## Assessing the improvement compared to the state-of-the-art (e.g., cytogenetic) risk:
            levels(riskLevel) <- c("Low", "Int", "High")
            formulaRc <- as.formula(paste("yTrain ~", paste(riskCol, "cond", sep=" + ")))
            r1 <- factor(trainDset[names(cond), riskCol], levels=riskLevel[levels(riskLevel)])
            c1 <- factor(cond, levels=c("Low", "Int", "High"))
            coxData1 <- as.data.frame(cbind(Risk=r1, cond=c1))
            fitCoxRc <- survival::coxph(formula=formulaRc, data=coxData1)
            capture.output("Coefficients of the Cox model using both Risk and cond:", 
                           file=outFile, append=TRUE)
            capture.output(summary(fitCoxRc), file=outFile, append=TRUE)
            textplot(summary(fitCoxRc)$coefficients)
            title("Cox coefficients")
        } ## if(doAddConfTable)
        dev.off()
        ##
        subsets[[subPowerSet]][["cond"]] <- cond

        ## #### Plot eigengenes for the selected modules:
        if(FALSE){ ## Not tested yet, TODO.
            Prognosis <- as.data.frame(trainDset[, eventCol])
            rownames(Prognosis) <- rownames(trainDset)
            colnames(Prognosis) <- eventCol
            Prognosis[Prognosis == 0] <- "DiseaseFree"
            Prognosis[Prognosis == 1] <- "Recurred/Progressed"
            par(mfrow=c(1, 1))
            png(file.path(resultPath, paste("eigen_", "-", subPowerSet, ".png", sep="")))
            pheatmap.type(Data=trainDset[, currentModule, drop=FALSE], annRow=Prognosis, 
                          cluster_cols=FALSE, show_rownames=FALSE)
            dev.off()
            capture.output(table(cond, trainDset[, eventCol]), 
                           file=outFile, append=TRUE)
        }
    }## End for.
    result[["subsets"]] <- subsets

    ## save pvalues
    result[["kmPVals"]] <- kmPVals

    ## Output:
    aft <- result
    save.if(aft, file=file.path(resultPath, "aft.RData"), verbose=verbose)
    return(result)
}
