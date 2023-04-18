\name{distanceToTss}
\alias{distanceToTss}

\title{Computes distance of loci to the closest TSS}

\description{
    Computes distance of the input loci to the closest Transcription Start Site
    (TSS). Also, creates a plot to identify the fraction of loci that have distance
    less than 1000 bp to TSS.}

\usage{
distanceToTss(usefuLoci, locus2oneGene, genesColName="Gene_Symbol", 
                        coordinatesColName="Genomic_Coordinate", 
                        lociColName="probeID", verbose=0)}

\arguments{
  \item{usefuLoci}{
    A list of character vectors. The name of each element of the list is a gene and
    the value is a vector of the corresponding loci, e.g, c("cg05116342", "cg15516558").}

  \item{locus2oneGene}{
    A data frame with at least four columns, namely: "probeID", "Gene_Symbol",
    "Chromosome", "Genomic_Coordinate". The rows are loci.}

  \item{genesColName}{
    A name of the column in the \code{locus2oneGene} dataframe that contains gene names.
    Default is set to "Gene_Symbol".}

  \item{coordinatesColName}{
    A name of the column in the \code{locus2oneGene} dataframe that contains genomic
    co-ordinate information. Default is set to "Genomic_Coordinate".}

  \item{lociColName}{
    A name of the column in the \code{locus2oneGene} dataframe that contains loci names.
    Default set to "probeID".}

  \item{verbose}{
    A numeric value that can set verbosity level. Default = 1 indicates minimal
    verbosity.}

}
\details{
In the iNETgrate package this function is called in the \code{\link{computeEigenloci}}
function. 
}

\value{
A list of:
  \item{distanceToClosesTssPerGene}{An integer vector named by the loci. If a
        locus is related to multiple genes, there would be multiple entries of the
        same loci, e.g., cg27665767.4, cg27665767.8, cg27665767.9, etc.}
  \item{distanceToClosesTss}{An integer vector named by the loci.}
  \item{transcriptLen}{A numeric vector of length of gene transcripts}
  \item{notFound}{A character vector of genes for which transcript length could
                  not be calculated.}
}

\seealso{
\code{\link{computeEigenloci}}
}

\examples{
data(toyCleanedAml)
genes1 <- c("FOXK1", "NTM", "MLPH", "SNTB1", "KCNA7", "XYLT1", "NF1")

## Creating locus2oneGene dataframe by calling createLocusGene function
createdLocusGene <- createLocusGene(locus2gene=toyCleanedAml$locus2gene,
                                     genesColName="Gene_Symbol",
                                     lociColName="probeID", verbose=0)
locus2oneGene <- createdLocusGene$locus2oneGene
subLocusGene <- locus2oneGene[is.element(locus2oneGene[ ,"Gene_Symbol"], genes1), ]

## A list of usefuLoci
usefuLoci <- split(subLocusGene$probeID, subLocusGene$Gene_Symbol)

computeDistance <- distanceToTss(usefuLoci=usefuLoci, 
                                           locus2oneGene=locus2oneGene, 
                                           genesColName="Gene_Symbol", 
                                           coordinatesColName="Genomic_Coordinate", 
                                           lociColName="probeID", verbose=1)
}

\keyword{manip}
