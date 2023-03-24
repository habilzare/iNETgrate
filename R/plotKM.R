plotKM <- function(inputTime, inputEvent, cond,
                   cond2Color=c("Low"='green',"Int"='blue',"High"='red'),
                   titleText=NULL, weight=rep(1, times=length(inputTime)),
                   doIntLow4pval=FALSE, plotFile=NULL, legendCex=1.7,
                   ylab1="Survival Probablity", pvalDigits=0, xmax1=Inf, xmin1=0, 
                   verbose=0, ...){
    ## Draws KM plots, 1 curve per each group of cases in the cond
    ## doIntLow4pval: If TRUE, intermediate cases are considered low-risk when computing the p-value.
    ##^Otherwise, they will be excluded.
    result <- list()
    message.if("Plotting KM...", verbose=verbose)
    result[["pVal"]] <- 1
    result[["survfitOutput"]] <- NULL
    doSkiPlot <- FALSE
    result[["doSkiPlot"]] <- doSkiPlot
    p1 <- NULL
    xmax2 <- xmax1
    if(!is.finite(xmax2)){ ##E.g., Inf
        xmax2 <- min(xmax1, max(inputTime, na.rm=TRUE))
    }
    xlim1 <- c(xmin1, xmax2)
    if(!is.null(plotFile)){
        png(plotFile, height=2*480, width=2*480, res=2.5*72)
        p1 <- plotKM(inputTime=inputTime, inputEvent=inputEvent, cond=cond, 
                     cond2Color=cond2Color,
                     titleText=NULL,  weight=weight,
                     doIntLow4pval=doIntLow4pval, plotFile=NULL,
                     legendCex=1, ylab1=ylab1, pvalDigits=pvalDigits,
                     xmax1=xmax1, xmin1=xmin1, verbose=verbose-2, ...)
        dev.off()
    }
    ## plotFile is NULL from here to the end of the function.

    ## QC:
    doAddPval <- !is.na(pvalDigits)
    inputEvent <- inputEvent[names(inputTime)]
    if(length(inputTime)<2){
        message.if("xlim is ignored because there is only 1 case.", verbose=verbose-2)
        xlim1 <- NULL
    }
    if(length(inputTime)!= length(inputEvent))
        stop("Do names match in inputTime and inputEvent?!")
    if(length(inputTime)==0){
        message.if("Length of inputTime is 0! KM plot is empty!", verbose=verbose-1)
        doSkiPlot <- TRUE
    }

    ## Making Data
    commonPatients <- intersect(names(inputTime), names(inputEvent))
    commonPatients <- intersect(commonPatients, names(cond))
    if(length(commonPatients)!=
       max(length(inputTime), length(inputEvent),length(cond))){
        message.if("inputTime, inputEvent, and cond must have the same names!",
                   verbose=verbose-1)
    }
    Data <- as.data.frame(cbind(Years=inputTime, Status=inputEvent))
    Data <- cbind(Data, cond=cond)
    if(nrow(Data)==0){
        message.if("Data has no rows! KM plot is empty!", verbose=verbose-1)
        doSkiPlot <- TRUE
    }

    
    if(verbose>1){
        print(table(cond))
    }
    if(doSkiPlot){
        message.if("Skipping this plot.", verbose=verbose)
        plot.new()
        return(result)
    }
    message.if("Surv...", verbose=verbose-3)
    yTrain <- Surv(time=inputTime, event=inputEvent, type="right")
    names(yTrain)<-names(inputTime)
    result[["yTrain"]] <- yTrain
    ##
    message.if("survfit...", verbose=verbose-3)
    survfitOutput <- survfit(formula=Surv(time=inputTime, event=inputEvent, type="right") ~ cond, data=Data)
    result[["survfitOutput"]] <- survfitOutput
    colorChart<-c()
    if(length(survfitOutput$strata)==0){
        if(any(!unique(cond) %in% names(cond2Color)))
            stop("cond has a value that does not exist in names of cond2Color!")
        colorChart[unique(cond)] <- cond2Color[unique(cond)]
    } else {
        for(survCurve in names(survfitOutput$strata)){
            survCurve<-gsub(survCurve,pattern=".*=",replacement="")
            colorChart[survCurve] <- cond2Color[survCurve]
        }
    }
    colorChartFitOrder <- colorChart
    ## Order: Low, Int, High
    colorChart <- colorChart[intersect(names(cond2Color),names(colorChart))]
    ## Add "Risk" for the legend:
    legendN <- paste0(names(colorChart), "-risk (n=")
    legendN <- paste0(legendN, table(cond)[names(colorChart)], ")")
    names(legendN) <- names(colorChart)

    ## Plot:
    if(!doSkiPlot){
        plotArgs <- list(x=survfitOutput, col=colorChartFitOrder, main=NULL, xlab="Time (year)",
                         ylab=ylab1, lwd=2, cex.axis=1.5, cex.lab=1.5, ...)
        if(!is.null(xlim1)){
            plotArgs <- c(plotArgs, list(xlim=xlim1))
        } else {
            plotArgs <- c(plotArgs, list(xmax=xmax1))
        }
        do.call(what="plot", args=plotArgs, quote=TRUE)
        legend("topright", col=colorChart, legendN[names(colorChart)], pch=19, cex=legendCex)
        title(titleText, cex.main= 1)
    } else {
        p1 <- plot.new()
    }
    
    
    ## Compute p-value for the model:
    if(length(intersect(cond, c("Low", "High")))==2){ ## Both Low and High are present,
        ##^they are used to compute p-value.
        message.if("Computing p-value...", verbose=verbose-4)## silent when verbose < 5.

        ## Log-Rank test to calculate p-value to put on legend of graph
        cond4Pval <- cond
        if(doIntLow4pval){
            cond4Pval[cond4Pval=="Int"] <- "Low"
        }
        ## We don't consider Medium (Int) risk in computing the p-value.
        survRes <- Surv(time=Data[cond4Pval!="Int","Years"],
                        event=Data[cond4Pval!="Int","Status"],type="right")
        survdiffRes <- survdiff(formula=survRes ~ cond4Pval[cond4Pval!="Int"])
        ##^ if the status and time in both conditions are the same, the above line
        ##^ leads to the following error;
        ## "Lapack routine dgesv: system is exactly singular: U[1,1]=0".
        pVal <- 1 - pchisq(survdiffRes$chisq, length(survdiffRes$n)-1)
        
        ## Try the following for better presentation:
        pVaLog <- round(log10(pVal), digits=pvalDigits)
        if(pVaLog < -2){
            pValToShow <- bquote("P-value:" ~ 10^.(pVaLog))
        } else {
            pValToShow <- bquote("P-value:" ~ .(round(pVal, digits=3)))
        }

        if(doAddPval){
            xPval <- 1.3+xlim1[1] ##1.05+xlim1[1]
            yPval <- 0.1*(xPval>=5) + (1.03)*(xPval<5)## when xPval is small, we put the pvalue at the top. 
            legend(x=xPval, y=yPval, box.lty=0, cex=legendCex*1.4, x.intersp=0, y.intersp=0,
                   legend=bquote(.(pValToShow)))
        }
        
        result[["pVal"]] <- pVal
    } ## if Both Low and High are present.
    return(result)
}
