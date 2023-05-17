\name{preprocessDnam}
\alias{preprocessDnam}

\title{Cleans DNA methylation data}

\description{
    Removes any loci with more than 50\%  missing beta values and
    imputes data for any locus that has less than 50 percent missing
    beta values (i.e., \code{NA}s are replaced with the mean of the
    corresponding locus).  Also, it identifies and removes
    loci that are identified as non-CpG, have Single Neucleotide Polymorphism
    (SNPs), or have missing chromosome information.}

\usage{
preprocessDnam(rawDnam, rawDnamSampleInfo=NULL, savePath, annLib="Auto", verbose=0)
}

\arguments{
  \item{rawDnam}{
    A matrix of beta values or a summarized experiment object similar to
    the output of \code{\link[TCGAbiolinks:GDCprepare]{GDCprepare}}.
    See 'Details' section for more.}

  \item{rawDnamSampleInfo}{
    A data frame that gives information about the DNA methylation samples with 
    at least two columns that contain patient IDs and corresponding sample IDs. 
    See 'Details' section for more.}

  \item{savePath}{
    A character vector that determines the path to the folder where plots can be saved.}

  \item{annLib}{
    A string that gives name of annotation (reference) library to be
    used. Default is "IlluminaHumanMethylation450kanno.ilmn12.hg19".}

  \item{verbose}{
    An integer that sets the level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\details{
If \code{rawDnam} is a matrix, the row names are probe IDs, column names are
patient IDs, and values are dnam beta values. In case \code{rawDnam} is of class 
\code{SummarizedExperiment}, it is assumed that sample information is a
part of \code{rawDnam} and hence, \code{rawDnamSampleInfo} can be set to
\code{NULL}. Otherwise \code{rawDnamSampleInfo} input is needed to get
an appropriate sampleInfo output.
}

\value{
A list with following components:
  \item{dnam}{A matrix of dnam beta values that is cleaned and imputed where 
              rownames are probe IDs and column names are patient IDs.}
  \item{locus2gene}{A dataframe of loci information  where rownames are probeIDs 
                    and at least 4 columns that contain information about probeIDs, 
                    geneIDs, genomic coordinates, and chromosome information.}
  \item{sampleInfo}{A dataframe of sample information where rownames are sample
                    IDs and at least two columns that contains patient
                    information and sample type. The rows are samples.}
}

\seealso{
\code{\link[TCGAbiolinks:GDCprepare]{GDCprepare}},
\code{\link[minfi:getAnnotation]{getAnnotation}}
}

\examples{
data(toyRawAml)
processeDnam <- preprocessDnam(rawDnam=toyRawAml$rawDnam,
                                rawDnamSampleInfo=NULL, savePath=tempdir(),
                                annLib="Auto", verbose=1)
}

\keyword{manip}
