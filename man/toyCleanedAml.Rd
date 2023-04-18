\name{toyCleanedAml}
\alias{toyCleanedAml}
\docType{data}
\title{A subset of cleaned TCGA-LAML data}

\description{
Survival data of 186 acute myeloid leukemia (AML) patients from The Cancer 
Genome Atlas (TCGA) managed by the NCI and NHGRI (Accession: phs000178.v11.p8)
is included here. Additionally, DNA methylation profile of 600 loci for 194 AML 
cases and gene expression profile of 396 genes (associated with loci) for 173 
AML cases are included here.}

\usage{data("toyCleanedAml")}

\format{
  A list}

\value{
  A list of \code{4} components.
  \item{genExpr}{A matrix of gene expression profile where rows and columns are 
                 genes Symbol, and sample IDs, respectively.}
  \item{locus2gene}{A dataframe where rows are probe IDs.}
  \item{dnam}{A matrix of DNA methylation profile where rows and columns are 
              probe IDs, and sample IDs, respectively.}
  \item{survival}{A dataframe where rows are patient IDs, and columns
                 provide survival information and vital status.}
}

\details{
    The original gene expression and DNA methylation data was produced using 
    Affymetrix U133 Plus 2 platform and Illumina Infinium HumanMethylation450 
    BeadChip platform, respectively. These TCGA--LAML data are aimed to serve 
    as a foundation for further investigations of AML pathogenesis, classification, 
    and risk stratification.}

\source{\url{https://portal.gdc.cancer.gov/projects/TCGA-LAML}}


\references{
    Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia
    (2013). New England Journal of Medicine 368.22: 2059-2074.}

\seealso{
  \code{\link{iNETgrate-package}},
  \code{\link{toyRawAml}}, \code{\link{toyComputEloci}}
}

\examples{


## This dataset can be regenerated using the following commands,
##which can take some time to run.
\donttest{
dataProject <- "TCGA-LAML"
dataPath <- file.path(tempdir(), "data", dataProject)
dir.create(dataPath, recursive=TRUE)
print(paste("Data will be saved in:", dataPath))

## Downloading data
rawData <- downloadData(dataProject=dataProject, savePath=dataPath)

## Cleaning data
print("Cleaning and preparing all data...")
riskCatCol <- "acute_myeloid_leukemia_calgb_cytogenetics_risk_category"
riskFactorCol <- "cytogenetic_abnormalities"

cleaned <- cleanAllData(genExpr=rawData$genExpr,
                         genExprSampleInfo=rawData$genExprSampleInfo, 
                         rawDnam=rawData$rawDnam, savePath=dataPath, 
                         annLib="Auto", clinical=rawData$clinical, 
                         riskCatCol=riskCatCol, riskFactorCol=riskFactorCol, 
                         riskHigh="Poor", riskLow="Favorable", otherLabel=NULL, 
                         verbose=1)

## Prepare toy data
print("Creating a sample of data!!")
toyData <- sampleData(Data=rawData, cleanData=cleaned, seed=1)

toyCleanedAml <- toyData$cleaned
}

data(toyCleanedAml)
class(toyCleanedAml)
}

\keyword{datasets}
