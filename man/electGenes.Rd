\name{electGenes}
\alias{electGenes}

\title{Selects input genes for iNETgrate analysis.}

\description{
    Selects a set of genes that have moderate or high correlation with
    survival time or vital status based on their expression or
    methylation levels. This is a  wrapper function that executes
    \code{\link{filterLowCor}} for filtering out low correlation genes
    followed by the \code{\link{computeUnion}} functions.}

\usage{
electGenes(genExpr, dnam, survival, savePath, event="Dead", locus2gene, 
            genesColName="Gene_Symbol", lociColName="probeID",
            ratio=c(genExpr=1/3, dnam=1), minCor=c(genExpr=0.2, dnam=0.2), 
            doAddTitle=TRUE, doAlLoci=FALSE, verbose=0)
}

\arguments{
  \item{genExpr}{
    A matrix of gene expression data where row names are gene IDs and
    column names are patient IDs or case IDs.}

  \item{dnam}{
    A matrix of beta values where row names are probe or loci IDs and column names 
    are patient IDs or case IDs.}

  \item{survival}{
    A matrix or data frame where rownames are patient IDs and must have minimum 
    of four columns namely: "Time" (to event), \code{event}, "Abnormality" and "Risk1".}

  \item{savePath}{
    A character vector determining the path to the folder or directory where the 
    plots and output are saved.}

  \item{event}{
    A string determining the name of the event column in the survival data. It stores 
    binary values 1 and 0 that correspond to the occurance or non-occurance of
    the event respectively.}

  \item{locus2gene}{
    A data frame or matrix with at least two columns that have information about
    probe IDs and the associated gene names. The rows are loci.}

  \item{genesColName}{
    A string determining the name of the column in \code{locus2gene} matrix that 
    carries gene ID information.}

  \item{lociColName}{
    A string determining the name of the column in \code{locus2gene} matrix that 
    carries probe ID information.}

  \item{ratio}{
    A named vector with two numeric values that determine the ratio of highly
    correlated genes and loci to be included in further
    analysis. The default is set to \code{1/3} for gene expression data
    and \code{1} for DNA methylation data.}

  \item{minCor}{
    A named vector with two numeric values that determine the minimum absolute 
    correlation cut--off value to filter data. Default set to \code{0.2} for 
    both gene expression and DNA methylation.}

  \item{doAddTitle}{
    A Boolean value, where \code{TRUE} (default) indicates appropriate 
    titles will be added to the plots.}

  \item{doAlLoci}{
    A Boolean value, where \code{FALSE} (default) indicates that union 
    gene set is composed of genes corresponding to selected (highly correlated) loci.}

  \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher 
    values produce more details of the computation.}
}

\value{
A list comprising of following objects:
  \item{selectedGenesExpr}{A character vector of gene IDs that show high
                           correlation with vital status or survival time
                           based on gene expression data.}
  \item{selectedLoci}{A character vector of loci or probe IDs that show high
                      correlation with vital status, survival time or both based
                      on DNA methylation data.}
  \item{selectedGenesDnam}{A character vector of gene IDs that have high
                           correlation with survival time or vital status 
                           based on DNA methylation data.}
  \item{unionGenes}{A character vector of gene IDs that have correlation with survival 
                    time or vital status based on expression levels or DNA
                    methylation data.}
  \item{unionLoci }{A character vector of loci or probe IDs that correspond to
                    the \code{unionGenes} set and also have appropriate methylation data.}
}

\seealso{
\code{\link{filterLowCor}}, \code{\link{computeUnion}}
}

\examples{

data(toyCleanedAml)

electedGenes  <- electGenes(genExpr=toyCleanedAml$genExpr, dnam=toyCleanedAml$dnam,
                             survival=toyCleanedAml$survival, savePath=tempdir(), 
                             locus2gene=toyCleanedAml$locus2gene, verbose=1)

}
\keyword{manip}
