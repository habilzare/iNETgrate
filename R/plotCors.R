plotCors <- function(corTimeData, corVitalData, corUnion, savePath,
                     corCutoff, whichData){
    if(length(whichData) !=1)
        stop("Please choose ONE correct input data type ",
             "(e.g., 'expression' or 'dnam')!!")
    xLabel <- "# of Loci"
    if(whichData == "expression")
        xLabel <- "# of genes"
    ## Plot cors:
    plotFile <- file.path(savePath, paste0(whichData, "_cors.png"))
    png(plotFile, height=500, width=3*500)
    par(mfrow=c(1,3))
    plot(y=sort(corTimeData[,"Time"]),x=seq_len(nrow(corTimeData)), 
         ylab="Correlation with survival time", xlab=xLabel, 
         cex.main=1.25, cex.lab=1.5, cex.axis=1.75, col="black")
    abline(h=corCutoff, col="orange", lty="dashed")
    abline(h=-corCutoff, col="orange", lty="dashed")
    abline(h=5, col="orange", lty="dashed")
    plot(y=sort(corVitalData[,"Vital"]),x=seq_len(nrow(corVitalData)), 
         ylab="Correlation with vital status", xlab=xLabel, 
         cex.main=1.25, cex.lab=1.5, cex.axis=1.75, col="black")
    abline(h=corCutoff, col="red", lty="dashed")
    abline(h=-corCutoff, col="red", lty="dashed")
    plot(corUnion, cex.main=3, cex.lab=1.5, cex.axis=1.75, col="black",
         main=paste("Correlation:", stats::cor(corUnion[,"Time"],corUnion[,"Vital"],
                                               use="complete.obs")))
    abline(h=corCutoff, col="red", lty="dashed")
    abline(h=-corCutoff, col="red", lty="dashed")
    abline(v=corCutoff, col="orange", lty="dashed")
    abline(v=-corCutoff, col="orange", lty="dashed")
    dev.off()
}
