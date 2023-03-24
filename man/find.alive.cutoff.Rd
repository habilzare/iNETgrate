\name{find.alive.cutoff}
\alias{find.alive.cutoff}

\title{Find Alive Cutoff}

\description{Finds the best (smallest) cutoff on hope where the ratio of alive cases, 
             for which time1 > until, is at least 'ratio'.}

\usage{
find.alive.cutoff(hope, time1, until, minRecall=0.2, risk="Low", 
                  Labels=NULL, doDebug=FALSE, resPath, verbose=0, 
                  doPlot=FALSE, lowLabel="Low", highLabel="High")}

\arguments{
  \item{hope}{
    A named vector where names are patient IDs and the values are the predicted
    survival time.}

  \item{time1}{
    A named vector where names are patient IDs and the values are survival time
    in years.The vector is of the same length as \code{hope}.}

  \item{until}{
    A numeric value e.g. 10 years for breast cancer. It will be ignored if 
    Labels is not \code{NULL}.}

  \item{minRecall}{
    Minimum recall value to classify patients based on their risk level.}

  \item{risk}{
    A string indicating the \code{risk} level for which cut--off is to be found. 
    Set to "Low" or "High".}

  \item{Labels}{
    A named vector where names are patient IDs same as time1 and hope and values
    are training risk groups: "High","Int", "Low". Alternatively, set 
    \code{until} to "left" to find the high--risk cut--off.}

  \item{doDebug}{
    A Boolean with a default value of \code{FALSE}.}

  \item{resPath}{
    A string describing the path to the folder to save results (plots). }

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}

  \item{doPlot}{
    A Boolean with a default value of \code{FALSE}.}

  \item{lowLabel}{
    A string describing the low risk--cases, e.g., "Low".}

  \item{highLabel}{
    A string describing the high--risk cases, e.g., "High".}
}

%% \details{}

\value{
A list which contains:
  \item{num}{A numeric value indicating the number of samples analysed.}
  \item{cutoff}{A number indicating the cutoff value for being high--risk or 
                low--risk based on survival time.}
  \item{highPrecision}{Precision for high--risk cut--off.}
  \item{lowPrecision}{Precision for low--risk cut--off.}
  \item{highRecall}{Recall for high--risk cut--off.}
  \item{lowRecall}{Recall for low--risk cut--off.}
}

\examples{
data(toyCleanedAml)
survival <- toyCleanedAml$survival
## Labels <- setNames(survival[, "Risk1"], nm=survival[, "PatientID"])
time1 <- setNames(survival[, "Time"]/365, nm=survival[, "PatientID"])
dummyTime <- time1
dummyTime[1:100] <- dummyTime[1:100]/2

aliveCutoff <- find.alive.cutoff(hope=dummyTime, time1=time1, until=1, 
                                 resPath=getwd(), doPlot=TRUE, verbose=2)

}

\keyword{survival}
