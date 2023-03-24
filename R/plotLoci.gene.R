plotLoci.gene <- function(distanceToTss, cutoff=1000, plotFile, verbose=0){
    ## Habil wrote this function to plot the genomic position of the selected loci with respect
    ## to the corresponding genes.
    ##distanceToTss: A vector of integers. Each element is the distance between a locus (prob) and
    ##the closest TSS of the nearby genes, e.g., from the output of the compute.distance.to.tss
    ##function, choose distanceToClosesTss.
    message.if(cat("The empirical cumulative distribution function", 
                "(CDF) to be plotted at:", 
                plotFile, sep="\n"), verbose=verbose-1)

    d1 <- abs(distanceToTss)
    d1 <- d1+1 ## log(0) is not defined.
    d1 <- d1[!is.na(d1)]
    xmax <- round(log10(max(d1, na.rm=TRUE)))+1
    ymax <- round(sum(d1<=cutoff)/length(d1), digits=2)
    if(!is.null(plotFile))
        png(plotFile, height=2*480, width=2*480, res=2.5*72)
    plot(ecdf(log10(d1)), main=NULL, xaxt="n",
         ylab="Cumulative probability", xlab="Distance to TSS (bps)")
    axis(1, at=0:xmax, labels=10^(0:xmax))
    mtext(side=1, at=log10(cutoff), col="red", text=cutoff)
    mtext(side=2, at=ymax, col="red", text=ymax)
    abline(v=log10(cutoff), col="red", lty=2)
    abline(h=ymax, col="red", lty=2)
    if(!is.null(plotFile))
        dev.off()

    return(d1-1)
}
