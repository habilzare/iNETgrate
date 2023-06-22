\name{filterLowCor}
\alias{filterLowCor}

\title{Filters out genes or loci}

\description{
    Filters genes or loci that have low correlation with survival time or vital
    status and selects highly correlated genes or loci for further analysis.}

\usage{
filterLowCor(Data, survival, whichData=c("expression", "dnam"), ratio=1/3, 
               minCor=0.2, event="Dead", savePath, verbose=0)}

\arguments{
  \item{Data}{
    A matrix of gene expression or beta values where row names are gene or probe
    IDs and column names are patient or case IDs.}

  \item{survival}{
    A matrix or data frame where rownames are patient IDs and must have minimum 
    of two columns namely "Time" (to event) and \code{event}.}

  \item{whichData}{
    A character value to describe input Data being filtered. E.g. expression for gene
    expression data.}

  \item{minCor}{
    Minimum absolute correlation cut--off value to filter data. Default set to
    \code{0.2}.}

  \item{ratio}{
    A ratio that determines highly correlated genes or loci to be included in further 
    analysis. Default set to \code{1/3}.}

  \item{event}{
    The name of the event column in the survival data. It stores binary values 
    \code{1} and \code{0} that correspond to the occurance or non-occurance of 
    the event, respectively.}

  \item{savePath}{
    Path to the directory where plots and output are to  be saved.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
A list comprising of following objects:
  \item{alList}{A character vector of complete list of gene or probe IDs from the
                input.}
  \item{selectedList}{A character vector of gene or probe IDs that show high
                      correlation with vital status or survival time.}
  \item{corValues}{A numeric vector of correlation values where names are
                   selected genes or loci.}
}

\seealso{
\code{\link{prepareSurvival}}
}

\examples{
data(toyCleanedAml)

## Gene expression
filteredGeneExpr <- filterLowCor(Data=toyCleanedAml$genExpr,
                                   survival=toyCleanedAml$survival, 
                                   savePath=tempdir(), whichData="expression")

## DNA methylation
filteredDnam <- filterLowCor(Data=toyCleanedAml$dnam,
                               survival=toyCleanedAml$survival, 
                               savePath=tempdir(), whichData="dnam", ratio=1)
}

\keyword{manip}
