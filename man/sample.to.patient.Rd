\name{sample.to.patient}
\alias{sample.to.patient}

\title{Maps sample IDs to corresponding patient IDs}

\description{
    Maps sample IDs in gene expression and DNA methylation data to corresponding
    patient IDs in clinical data and renames columns of the matrix with patient IDs. 
    The function also provides options to exclude or include FFPE samples, undesired 
    sample types, or duplicate samples.}

\usage{
sample.to.patient(sampleInfo, sampleCol="barcode", patientCol="patient", 
                  sampleTypeCol="shortLetterCode", sampleTypesIn=NULL, isFfpe=NULL, 
                  doRemoveDup=TRUE, verbose=0)}

\arguments{
  \item{sampleInfo}{
    A data frame or matrix that where rows are the samples with at least 
    two columns that determines patient IDs and corresponding sample IDs.}

  \item{sampleCol}{
    A string that determines the name of the column that contains sample IDs
    in \code{sampleInfo} data frame.}

  \item{patientCol}{
    A string that determines the name of the column that contains patient IDs
    in \code{sampleInfo} data frame.}

  \item{sampleTypeCol}{
    A string that determines the name of the column that contains sample type
    information in \code{sampleInfo} data frame e.g., "shortLetterCode" in
    the \code{toyRawAml} data. Set it to NULL to disable filtering based
    on sample types. See \code{sampleTypesIn}.}

  \item{sampleTypesIn}{
    A character vector that identifies the sample types to be included for
    iNETgrate analysis e.g., c("TP", "NT") are sample types in TCGA data where
    TP and NT indicate Primary Tumor and Non-Tumor adjacent samples respectively.}
 
  \item{isFfpe}{
    A character vector that takes Boolean values to indicate inclusion or
    exclusion of FFPE sample e.g., \code{c(FALSE, TRUE)} will include
    all samples, \code{TRUE} will include only FFPE samples, and
    \code{FALSE} will exclude FFPE samples.
    Set it to \code{NULL} to disable FFPE sample filtering.}

  \item{doRemoveDup}{
    A Boolean value that indicates inclusion or exclusion of duplicate samples.
    E.g., \code{TRUE} (default) removes duplicate samples using
    lexicographic sorting of the sample IDs as described in
    \code{\link{find.tcga.duplicates}}}. 

   \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\details{
    Frozen samples are preferred over FFPE samples for sequencing, and hence,
    removal of FFPE data is preferred for any expression or methylation
    analysis.}

\value{
A list with following components:
  \item{sample2patient}{A named vector where names are sample IDs and values are
                        patient IDs.}
  \item{patient2type}{A named vector where names are tagged patient IDs and
                      values are the sample type e.g., "TB" is the sample type
                      in the \code{toyRawAml} data}
  \item{patient2type}{A list of excluded sample IDs}
}

\seealso{
\code{\link{clean.allData}}, \code{\link{find.tcga.duplicates}} 
}

\examples{
data(toyRawAml)

s2pExpr <- sample.to.patient(sampleInfo=toyRawAml$genExprSampleInfo, 
                             sampleCol="barcode", patientCol="patient", 
                             sampleTypeCol="shortLetterCode",
                             sampleTypesIn=NULL, isFfpe=NULL, 
                             doRemoveDup=TRUE, verbose=1)
}
\keyword{manip}
