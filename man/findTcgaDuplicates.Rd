\name{findTcgaDuplicates}
\alias{findTcgaDuplicates}

\title{Identify and remove duplicate samples}

\description{
    Identifies sample IDs in gene expression and DNA methylation data that belong to
    same patients or obtained from same cell types or have been obtained or
    processed using same experimental techniques. The samples are removed using
    lexicographic sorting of the sample IDs as generally used in sorting TCGA
    samples.}

\usage{
findTcgaDuplicates(sampleInfo, sampleCol="barcode", patientCol="patient", 
                     sampleTypeCol="shortLetterCode", doTag=FALSE, verbose=0)
}

\arguments{
  \item{sampleInfo}{
    A data frame or a matrix that gives information about the samples with at least 
    two columns that contain patient IDs and corresponding sample IDs. The rows
    are samples.}

  \item{sampleCol}{
    A string that determines the name of the column in \code{sampleInfo} data frame 
    that contains sample IDs.}

  \item{patientCol}{
    A string that determines the name of the column in \code{sampleInfo} data frame 
    that contains patient IDs.}

  \item{sampleTypeCol}{
    A string that determines the name of the column in \code{sampleInfo} data frame 
    that contains sample type information e.g., "shortLetterCode" in toyRawAml data. 
    Set it to \code{NULL} to disable filtering based on sample types.}

  \item{doTag}{
    A Boolean value that indicates if patient IDs in patientCol will be tagged
    with sample type from sampleTypeCol e.g., \code{TRUE} tags patient IDs with sample
    types that prevents removing duplicate samples from same patients with
    different sampling source and vice versa.}
    
   \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
A list with following components:
  \item{sampleInfo}{A data frame or a matrix with the same number of columns as
                    the input \code{sampleInfo}, but fewer rows. 
                    Rows of duplicate samples are removed.}
  \item{keptDups}{A character vector of duplicate sample IDs that are kept for further
                  analysis.}
  \item{excludeDups}{A character vector of duplicate sample IDs that are
                     excluded from further analysis.}
}

\seealso{
\code{\link{sample2pat}}
}

\examples{
data(toyRawAml)

removeDups <- findTcgaDuplicates(sampleInfo=toyRawAml$genExprSampleInfo, 
                                   sampleCol="barcode", patientCol="patient", 
                                   sampleTypeCol="shortLetterCode",
                                   doTag=TRUE, verbose=1)
}
\keyword{manip}
