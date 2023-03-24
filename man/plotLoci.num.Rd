\name{plotLoci.num}
\alias{plotLoci.num}

\title{Plot number of loci per gene.}

\description{Creates plots of the number of loci per gene.}

\usage{
plotLoci.num(locus2gene, genesColName="Gene_Symbol", lociColName="probeID", 
             selectedLoci="Auto", plotFile, doAddTitle=TRUE, 
             xlab1="Number of loci per gene", ylab1="Cumulative probability",
             verbose=0)
}

\arguments{
  \item{locus2gene}{
    A matrix or a dataframe with at least two columns that contain information
    about probe ID (e.g., cg00000029), and corresponding gene symbol (e.g.,
    TSPY4). The rows are loci.}

  \item{genesColName}{
    A string that determines the gene IDs column in locus2gene.}

  \item{lociColName}{
    A string that determines the probe IDs column in locus2gene.}

  \item{selectedLoci}{
    A character vector of probe IDs that are selected to be plotted.}

  \item{plotFile}{
    A string that determines the file path to save plots.}

  \item{doAddTitle}{
    Accepts a Boolean, where \code{TRUE} (default) indicates appropriate titles will be
    added to the plots.}

  \item{xlab1}{
    A string describing the X--axis label for the plots.}

  \item{ylab1}{
    A string describing the Y--axis label for the plots.}

  \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
It is a list with following components:
  \item{noGenes}{A character vector of probe IDs that have no associated gene
                 assignement.}
  \item{l2g}{A sorted table of number of loci associated with each gene.}
}

\examples{
data(toyCleanedAml)

savePath <- file.path(getwd(), "plotLociPerGene.png") 

plottedLoci <- plotLoci.num(locus2gene=toyCleanedAml$locus2gene, genesColName="Gene_Symbol", 
                            lociColName="probeID", selectedLoci="Auto", plotFile=savePath, 
                            doAddTitle=TRUE, xlab1="Number of loci per gene", 
                            ylab1="Cumulative probability")
}

\keyword{dplot}% use one of  RShowDoc("KEYWORDS")
