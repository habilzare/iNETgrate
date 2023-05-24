\name{toyRawAml}
\alias{toyRawAml}
\docType{data}
\title{A subset of TCGA-LAML data}

\description{
Clinical data of 200 acute myeloid leukemia (AML) patients from The Cancer 
Genome Atlas (TCGA) managed by the NCI and NHGRI (dbGaP Accession: phs000178.v11.p8)
is included here. Additionally, DNA methylation profile of 1800 loci for 194 AML 
cases and gene expression profile of 1077 genes (associated with loci) for 173 
AML cases are included here.}

\usage{data("toyRawAml")}

\format{
  A list}

\value{
  It is a list of \code{4} components.
  \item{genExpr}{A matrix of gene expression profile where rows and columns are 
                 genes Symbol, and sample IDs, respectively.}
  \item{genExprSampleInfo}{A dataframe where rows are sample IDs, and columns
                 provide sample information.}
  \item{rawDnam}{A matrix of DNA methylation profile where rows and columns are 
                 probe IDs, and sample IDs, respectively.}
  \item{clinical}{A dataframe where rows are patient IDs, and columns
                 provide clinical information.}
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
  \code{\link{toyCleanedAml}}, \code{\link{toyComputEloci}}
}

\examples{

## This dataset can be regenerated using the following commands,
##which can take some time to run. 
\donttest{
dataProject <- "TCGA-LAML"
dataPath <- file.path(tempdir(), "data", dataProject)
dir.create(savePath)

## Downloading data
rawData <- downloaData(dataProject=dataProject, savePath=dataPath)

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
print("Creating sample of data...")
set.seed(seed=1)
toyData <- sampleData(Data=rawData, cleanData=cleaned, seed=1)
toyRawAml <- toyData$rawData
}

data(toyRawAml)
class(toyRawAml)
}

\keyword{datasets}
