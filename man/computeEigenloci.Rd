\name{computeEigenloci}
\alias{computeEigenloci}

\title{Computes an eigenloci for each gene}

\description{
    An eigenloci of a gene is a weighted average (PC1) of beta values of
    the most correlated loci (i.e., core) of that gene.}
  
\usage{
computeEigenloci(dnam, locus2gene, geNames, Labels, Label1, 
                  Label2, plotPath=NULL, genesColName="Gene_Symbol",
                  coordinatesColName="Genomic_Coordinate", lociColName="probeID", 
                  dnamGene=NULL, doDebug=FALSE, verbose=0)}

\arguments{
  \item{dnam}{
    A matrix of beta values where row names are probe IDs and column names are patient IDs.}

  \item{locus2gene}{
    A data frame with loci information where row names are probe IDs and must have
    at least 4 columns that contain probe IDs, gene IDs, genomic coordinates, and 
    chromosome information.}

  \item{geNames}{
    A character vector of selected genes for which eigenloci needs to be obtained
    e.g. "unionGenes", the output of \code{\link{electGenes}}. If
    \code{NULL} (default), eigenloci for all the genes in the gene IDs
    column of the \code{locus2gene} dataframe will be calculated.}

  \item{Labels}{
    A named character vector where names are tagged patient IDs and values are 
    group label. It is needed for sample balancing when computing the
    eigenloci. See \code{Label1} and  \code{Label2}.}

  \item{Label1}{
    One of the values from \code{Labels} vector that defines a patient group to be used to
    balance data e.g., low--risk group defined as "Favorable" in \code{toyAmlData}.}

 \item{Label2}{
    One of the values from \code{Labels} vector that defines a patient group to be used to
    balance data e.g., high--risk group defined as "Poor" in \code{toyAmlData}.}

  \item{plotPath}{
    A string that determines the path to save plots.}

  \item{genesColName}{
    A string determining the name of the column in the \code{locus2gene} data frame 
    that stores gene IDs.}

  \item{coordinatesColName}{
    A string determining the name of the column in the \code{locus2gene} data frame
    that stores genomic co-ordinates.}

  \item{lociColName}{
     A string determining the name of the column in the \code{locus2gene} data frame
     that stores probeIDs.}

  \item{dnamGene}{
    A character vector of selected genes based on correlating survival time with
    dnam data.}

  \item{doDebug}{
    A Boolean that aids in code debugging by experts only. Otherwise,
    leave it as \code{FALSE} (default).}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\details{
  For each genes, first its core (i.e., the most connected cluster of its
  loci) is identified using the \code{\link{findCore}} function. Then, an eigenloci
  for that gene is computed using PCA on the core loci, which optimizes the weights
  of each locus in the core.
}

\value{
  \item{eigenloci}{A matrix of eigenloci values where row names are
                   patient IDs and column names are gene IDs.}
  \item{usefuLoci}{A list of character vectors, where a gene is an item of the
                   list that stores the corresponding loci information.}     
  \item{lociPigen}{A list of Pigengene objects. Names are genes.}
  \item{distanceToTss}{A list of distances of all selected loci to 
                       Transcription Start Sites (TSS).}
  \item{distanceToTssDnam}{A list of distances of loci selected based on dnam data to TSS.}
}

\references{
Amir Foroushani \emph{et al.} (2016) Large-scale gene network analysis reveals the significance
of extracellular matrix pathway and homeobox genes in acute myeloid leukemia: an
introduction to the Pigengene package and its applications, \emph{BMC MedicalGenomics}.
}

\seealso{
\code{\link[Pigengene:compute.pigengene]{compute.pigengene}}, 
\code{\link[Pigengene:project.eigen]{project.eigen}},
\code{\link{distanceToTss}}, \code{\link{electGenes}}
}

\examples{
data(toyCleanedAml)
plotPath <- file.path(tempdir(), "plots")
print(paste("Results will be saved in:", plotPath))
dir.create(plotPath, showWarnings=FALSE)
unionGenes <- c("BAI1", "CAPNS1", "CHD1", "CYC1", "EXT1", "FYB", "GFRA1", "HLX", 
                "IDH3B", "PSMD8", "SGTA", "SLC1A5", "SSFA2", "TFAM", "TGIF1", 
                "UBC", "ZNF70", "ST8SIA4", "ZNF485", "FOXN3", "HCG4", "FOXK1", 
                "C2orf44")

survivalAml <- toyCleanedAml$survival
patientLabel <- setNames(as.character(survivalAml[,"Risk1"]), rownames(survivalAml))
inBoth <- intersect(colnames(toyCleanedAml$dnam), names(patientLabel))
dnamData <- toyCleanedAml$dnam[ ,inBoth]
eigenLoci <- computeEigenloci(dnam=dnamData, locus2gene=toyCleanedAml$locus2gene,
                               geNames=unionGenes,
                               Labels=patientLabel[is.element(names(patientLabel), inBoth)], 
                               Label1="Low", Label2="High", plotPath=plotPath,
                               verbose=1) 
}

\keyword{multivariate}
