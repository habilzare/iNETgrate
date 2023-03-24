\name{plotLoci.gene}
\alias{plotLoci.gene}

\title{Plots distance of loci to genes}

\description{Generates plots describing distance of a loci from a Transcription 
            Start Site (TSS).}
\usage{
plotLoci.gene(distanceToTss, cutoff=1000, plotFile, verbose=0)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distanceToTss}{
    A numeric vector, where each element is the distance between a locus (probe) and
    the closest Transcription Start Site (TSS) of the nearby genes.
    I.e., from the full output of the \code{\link{compute.distance.to.tss}} function, 
    only \code{distanceToClosesTss} was selected.}

  \item{cutoff}{
    A numeric value that determines the maximum distance (in basepairs) between
    probe and TSS.}

  \item{plotFile}{
    A character string that gives path to the folder where plot will be saved.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
The distance matrix. Also, a plot is saved at plotFile.}

\examples{
distanceToTss <- c(2460, 1498, 17182, 7807, 237, 56815, 368948)
names(distanceToTss) <- c("cg25539759", "cg25408497", "cg20716080", "cg07917842", 
                          "cg07312445", "cg04774620", "cg06129210")
savePath <- file.path(getwd(), "plotDistanceToTss.png") 

plotted <- plotLoci.gene(distanceToTss=distanceToTss, cutoff=1000,
                         plotFile=savePath, verbose=0)
}

\keyword{dplot}% use one of  RShowDoc("KEYWORDS")
