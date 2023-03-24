\name{toyComputEloci}
\alias{toyComputEloci}
\docType{data}
\title{A subset of weighted DNA methylation profile of TCGA-LAML data}

\description{
    Weighted average of DNA methylation profile of 180  acute myeloid leukemia (AML) cases 
    from TCGA-LAML data.}

\usage{data("toyComputEloci")}

\format{
  A list}

\value{
  A list of \code{4} components.
  \item{eigenloci}{A matrix of eigenloci values where row names are
                   patient IDs, and column names are gene IDs.}
  \item{usefuLoci}{A list of character vectors, where a gene is an item of the
                   list that stores the corresponding loci information.}     
  \item{lociPigen}{A list of Pigengene objects. Names are genes.}
  \item{distanceToTss}{A list of distances of all selected loci to 
                       Transcription Start Site(TSS).}
}

\source{\url{https://portal.gdc.cancer.gov/projects/TCGA-LAML}}


\references{
    Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia
    (2013). New England Journal of Medicine 368.22: 2059-2074.}

\seealso{
  \code{\link{iNETgrate-package}},
  \code{\link{toyCleanedAml}}, \code{\link{toyEigengenes}}
}

\examples{
data(toyComputEloci)

## How to reproduce the data:
data(toyCleanedAml)
survivalAml <- toyCleanedAml$survival
plotPath <- file.path(getwd(), "plots")
dir.create(plotPath, showWarnings=FALSE)

print("Electing genes...")
elected <- elect.genes(genExpr=toyCleanedAml$genExpr, dnam=toyCleanedAml$dnam,
                       survival=survivalAml, savePath=getwd(), 
                       locus2gene=toyCleanedAml$locus2gene, verbose=1)

print("Computing eigenloci!!")
patientLabel <- setNames(as.character(survivalAml[,"Risk1"]), nm=rownames(survivalAml))
inBoth <- intersect(colnames(toyCleanedAml$dnam), names(patientLabel))
computedEloci <- compute.eigenloci(dnam=toyCleanedAml$dnam[ ,inBoth], 
                                   locus2gene=toyCleanedAml$locus2gene, 
                                   geNames=elected$unionGenes,
                                   Labels=patientLabel[is.element(names(patientLabel), inBoth)], 
                                   plotPath=plotPath, Label1="High", Label2="Low", 
                                   verbose=1)

class(computedEloci)
}

\keyword{datasets}
