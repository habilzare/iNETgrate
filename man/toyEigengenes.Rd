\name{toyEigengenes}
\alias{toyEigengenes}
\docType{data}
\title{A subset of weighted average of gene expression present in each module}

\description{
    Weighted average of gene expression of 158  acute myeloid leukemia (AML) cases 
    from TCGA-LAML data that is distributed into 4 feature sets (modules).}

\usage{data("toyEigengenes")}

\format{
  A matrix}

\value{
  It is a \code{158*4} numeric matrix.}

\details{
    The columns and rows are named according to the feature names (modules),
    and patient IDs, respectively.}

\source{\url{https://portal.gdc.cancer.gov/projects/TCGA-LAML}}


\references{
    Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia
    (2013). New England Journal of Medicine 368.22: 2059-2074.}

\seealso{
  \code{\link{iNETgrate-package}}, \code{\link{computeEigengenes}},
  \code{\link{toyCleanedAml}}, \code{\link{toyComputEloci}}
}

\examples{
data(toyEigengenes)

## How to reproduce data:
data(toyCleanedAml)
data(toyComputEloci)
eigenloci <- toyComputEloci$eigenloci 
patientLabel <- as.matrix(toyCleanedAml$survival)[,"Risk1"]

# Path to the network
netPath <- file.path(tempdir(), "net")
dir.create(netPath, showWarnings=FALSE)
madeNetwork <- makeNetwork(genExpr=toyCleanedAml$genExpr, eigenloci=eigenloci,
                            geNames=colnames(eigenloci), mus=c(0.6), 
                            outPath=netPath, verbose=0)

eigengenes <- computeEigengenes(genExpr=toyCleanedAml$genExpr, 
                                  eigenloci=eigenloci, netPath=netPath, 
                                  geNames=colnames(eigenloci), Labels=patientLabel, 
                                  Label1="Low", Label2="High", 
                                  mus=c(0.6), survival=toyCleanedAml$survival,
                                  mu2modules=madeNetwork$mu2modules)

class(eigengenes)
}

\keyword{datasets}
