\name{sampleData}
\alias{sampleData}

\title{Sample Data}

\description{
    Samples data for quick running of the pipeline.}

\usage{
sampleData(Data, cleanData, numbeRow=300)}

\arguments{
  \item{Data}{
    A list of input data same as the output of \code{downloaData} function.}

  \item{cleanData}{
    A list of cleaned data same as the output of \code{\link{cleanAllData}} function.}

  \item{numbeRow}{
    An integer that indicates the number of rows to be sampled.}

}

\value{
A list with following elements:
  \item{rawData}{A list of sampled gene expression matrix, gene expression
                 sample information, rawDnam data, and clinical data.}
  \item{cleaned}{A list of sampled cleaned gene expression matrix, locus2gene
                 dataframe, cleaned Dnam data, and survival data.}
}

\note{
  Use \code{\link[base:set.seed]{set.seed}} before calling this function so that you
  will be able to reproduce the same sample in the future if needed.
}

\seealso{
  \code{downloaData}, \code{\link{cleanAllData}}
}

\examples{
data(toyRawAml)
data(toyCleanedAml)
sampleData(Data=toyRawAml, cleanData=toyCleanedAml, numbeRow=100)
}

\keyword{manip}
