\name{clean.allData}
\alias{clean.allData}

\title{Cleans and preprocesses data}
\description{
    This is a wrapper function that executes
    \code{\link{preprocess.dnam}} and
    \code{\link{prepare.survival}} functions to preprocess and impute DNA methylation 
    and clinical data respectively. It also creates the
    \code{\link{sample.to.patient}} mappings for gene expression and DNA
    methylation data.}

\usage{
clean.allData(genExpr, genExprSampleInfo, rawDnam, rawDnamSampleInfo=NULL, 
              sampleCol="barcode", patientCol="patient", sampleTypeCol="shortLetterCode", 
              sampleTypesIn=NULL, savePath, isFfpe=c(FALSE), annLib="Auto",
              clinical, patientIDCol="bcr_patient_barcode", eventCol="vital_status", 
              event="Dead", timeCol="days_to_last_followup", riskFactorCol,
              riskCatCol=NULL, riskHigh=NULL, riskLow=NULL, otherLabel=NULL, verbose=0)
}

\arguments{
  \item{genExpr}{
    A matrix of gene expression data where rows are genes and columns are samples.}

  \item{genExprSampleInfo}{
    A data frame that gives information about the gene expression samples. It
    contains at least two columns determining patient IDs and the corresponding 
    sample IDs. The row names may or may not be patient IDs.}

  \item{rawDnam}{
    A matrix of beta values or a summarized experiment object similar to
    the output of \code{\link[TCGAbiolinks:GDCprepare]{GDCprepare}}.
    Rows are genes and columns are samples.}

  \item{rawDnamSampleInfo}{
    A data frame that gives information about the DNA methylation samples. It
    contains at least two columns determining patient IDs and the corresponding 
    sample IDs. The row names may or may not be patient IDs.}

  \item{sampleCol}{
    The name of the column in the sample information data frames that contains sample IDs.}

  \item{patientCol}{
    The name of the column in the sample information data frames that contains patient IDs.}

  \item{sampleTypeCol}{
    The name of the column in the sample information data frames that
    contains sample type information. See \code{sampleTypesIn}.}

  \item{sampleTypesIn}{
    A character vector determining the sample types to be included in the
    iNETgrate analysis e.g., \code{c("TP", "NT")}  in TCGA data, where TP and
    NT indicate Primary Tumor and Non-Tumor adjacent samples,
    respectively.}
 
  \item{savePath}{
    A string determining the path to the folder to save plots.}

  \item{isFfpe}{
    A character vector that takes Boolean values to indicate inclusion or
    exclusion of FFPE sample. \code{c(FALSE, TRUE)} will include all
    samples, TRUE will include only FFPE samples, and FALSE will exclude
    FFPE samples. Set it to \code{NULL} to disable it.}

  \item{annLib}{
    A string that determines annotation (reference) library to be
    used. Default is "IlluminaHumanMethylation450kanno.ilmn12.hg19".}

  \item{clinical}{
    A matrix or a data frame of clinical data, with minimum four columns that gives 
    information about patient IDs, event occurance, survival or follow up time, 
    and risk factors. The rownames are patient IDs.}

  \item{patientIDCol}{
    A string determining the column in clinical data that lists patient
    or case IDs, e.g., "bcr_patient_barcode" in sample AML data}

  \item{eventCol}{
    A string determining the column in clinical data that lists the
    event. The event column must be a binary column e.g., "vital_status"
    column in sample AML data where two possible events are "Dead"
    and "Alive".}

  \item{event}{
    A string determining an event of interest from \code{eventCol} in clinical data. 
    E.g., "Dead" is the event of interest in the example data.}

  \item{timeCol}{
    A string determining the column in clinical data storing survival time
    information, e.g. "days_to_last_followup" in sample AML data}

  \item{riskFactorCol}{
    A string determining the column in clinical data that describes the risk factors
    associated to risk levels.}

  \item{riskCatCol}{
    A string determining the column in clinical data that lists risk levels of the disease.}

  \item{riskHigh}{
    A character vector that describes high risk group of patients in clinical data. 
    See section "Details".}

  \item{riskLow}{
    A character vector that describes low risk group of patients in clinical data. 
    See section "Details".}

  \item{otherLabel}{
    A vector of other risk levels where names are the same as row names of clinical data.}

  \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\details{
    When \code{riskCatCol} is \code{NULL}, \code{riskHigh} and \code{riskLow} 
    are expected to have values that belong to \code{riskFactorCol}. When 
    \code{riskCatCol} is defined, \code{riskHigh} and \code{riskLow} are expected 
    to have values that belong to \code{riskCatCol}. All \code{riskCatCol}, 
    \code{riskHigh} and \code{riskLow} cannot be \code{NULL} at the same
    time.}

\value{
A list with following components:
  \item{survival}{A data frame with a minimum of five columns namely: "PatientID", 
                  \code{event}, "Time", "Risk1", and "Abnormality". 
                  The row names are the same as "PatientID".}
  \item{genExpr}{A data frame of gene expression data where rows are genes and 
                 columns are (tagged) patient IDs.}
  \item{sample2patient}{A list of named vectors.}
  \item{patient2type}{A list of named vectors.}
  \item{dnam}{A matrix of dnam beta values that is cleaned and imputed where 
              rownames are probe IDs and column names are (tagged) patient IDs.}
  \item{locus2gene}{A data frame where rownames are probeIDs and at least 4 
                    columns that contain information about probeIDs, geneIDs, 
                    genomic coordinates, and chromosome information.}
}

\seealso{
\code{\link{prepare.survival}}, \code{\link{preprocess.dnam}}, 
\code{\link{sample.to.patient}}
}

\examples{
data(toyRawAml)
riskCatCol <- "acute_myeloid_leukemia_calgb_cytogenetics_risk_category"
riskFactorCol <- "cytogenetic_abnormalities"
cleaned <- clean.allData(genExpr=toyRawAml$genExpr,
                         genExprSampleInfo=toyRawAml$genExprSampleInfo, 
                         rawDnam=toyRawAml$rawDnam, savePath=getwd(), 
                         clinical=toyRawAml$clinical, riskCatCol=riskCatCol, 
                         riskFactorCol=riskFactorCol, riskHigh="Favorable", 
                         riskLow="Poor", verbose=1)
}

\keyword{manip}
