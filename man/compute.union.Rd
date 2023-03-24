\name{compute.union}
\alias{compute.union}

\title{Computes union gene set}

\description{
    Makes a union set of genes that have high correlation with survival time or
    vital status based on their expression or methylation levels.}

\usage{
compute.union(Genes, selectedGenes, loci, selectedLoci, locus2gene,
              doAlLoci=FALSE, genesColName="Gene_Symbol", lociColName="probeID", 
              verbose=0)}

\arguments{
  \item{Genes}{
    A character vector of all gene IDs from the gene expression data.}

  \item{selectedGenes}{
    A character vector of selected gene IDs that show high correlation with
    survival time or vital status based on their expression values.}

  \item{loci}{
    A character vector of all loci or probe IDs obtained after preprocessing
    methylation as detailed in \code{\link{preprocess.dnam}}.}

  \item{selectedLoci}{
    A character vector of loci or probe IDs that have high correlation with
    survival time or vital status based on their beta values.}

  \item{locus2gene}{
    A data frame or a matrix with at least two columns that contain information
    about probe IDs and the corresponding gene IDs. The rownames are probe IDs.}

  \item{doAlLoci}{
    A Boolean value with a default value of \code{FALSE}. If
    \code{TRUE}, a gene with any levels of DNA methylation loci
    will be included in the output union even if they do not correlate with
    survival data.}

  \item{genesColName}{
    A string determining the name of the column in \code{locus2gene} matrix that 
    carries gene ID information.}

  \item{lociColName}{
    A string determining the name of the column in \code{locus2gene} matrix that 
    carries probe ID information.}

  \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher 
    values produce more details of the computation.}
}

\value{
A list comprising of following objects:
  \item{unionGenes}{A character vector of selected genes that show high correlation
                    with survival time or vital status based on expression or
                    methylation data.}
  \item{unionLoci}{A character vector of loci that are associated with
                   unionGenes set and also have appropriate methylation data.}
  \item{dnamSelectGenes}{A character vector of selected genes that have 
                         correlation with survival time or vital status based on 
                         mehtylation data only.}
}

\seealso{
\code{\link{create.locusGene}}, \code{\link{elect.genes}}, 
\code{\link{preprocess.dnam}}
}

\examples{
data(toyCleanedAml)

highCorGenes <- c("TRPV4", "CHN2", "WFDC1", "XPO1", "HTR7", "ABCC1", "TMCO3", 
                  "FBXO22", "VWCE", "KLHL36", "APBA2", "IRF2BP2")
highCorLoci <- c("cg00044050", "cg00086266", "cg00112238", "cg00321478", 
                 "cg00371239", "cg00375235", "cg00390484", "cg00424690",
                 "cg00448143", "cg00452755")
unionSet <- compute.union(Genes=rownames(toyCleanedAml$genExpr),
                          loci=rownames(toyCleanedAml$dnam),
                          selectedGenes=highCorGenes, selectedLoci=highCorLoci,
                          locus2gene=toyCleanedAml$locus2gene, verbose=2)
}
\keyword{manip}
