\name{plotKM}
\alias{plotKM}

\title{Survival (Kaplan--Meier estimator) plots}

\description{
    Plots KM (survival) plots, one curve per each risk group cases.}

\usage{
plotKM(inputTime, inputEvent, cond, 
       cond2Color=c("Low"='green', "Int"='blue', "High"='red'), 
       titleText=NULL, weight=rep(1, times=length(inputTime)),  
       doIntLow4pval=FALSE, plotFile=NULL, legendCex=1.7, 
       ylab1="Survival Probablity", pvalDigits=0, xmax1=Inf, xmin1=0, 
       verbose=0, ...)}

\arguments{
  \item{inputTime}{
    A named numeric vector where names are patient IDs and values are survival 
    time (in years).}

  \item{inputEvent}{
    A named binary vector where names are patient IDs and values are \code{0}
    and \code{1}. Here, \code{1} indicates the occurance of the event and vice versa.}

  \item{cond}{
    A named character vector, where names are patient IDs and values are their
    risk category the patientID belongs to.}

  \item{cond2Color}{
    A named character vector, that assigns color to risk groups in KM--plots. Names 
    must be in the cond vector, and values are colors. 
    e.g, cond2Color=c("Low"='green', "Int"='blue', "High"='red')}

  \item{titleText}{
    A character string that assigns title of the plot.}

  \item{weight}{
    A numeric vector, same as length of inputTime vector.}

  \item{doIntLow4pval}{
    Accepts a Boolean value, where \code{FALSE} (default) indicates intermediate cases 
    are excluded when computing the p--value. If \code{TRUE}, intermediate cases are
    considered as low--risk cases for p--value computation.}

  \item{plotFile}{
    A character string that determines path to save the plots.}

  \item{legendCex}{
    A numeric value that determines graphical location of legend in the plots.}

  \item{ylab1}{
    A character string that determines the Y--axis label in the plot.}

  \item{pvalDigits}{
    A numeric value that gives number of decimal digits to be printed on the K--M
    plots for p--value.} 

  \item{xmax1}{
    A numeric value giving the maximum time (in years) for plotting the K--M curves.}

  \item{xmin1}{
    A numeric value giving the minimum time (in years) for plotting the K--M curves.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}

  \item{...}{
    User can pass any more arguements for creating better plots.}
}

\value{
A list with following components:
  \item{pVal}{A numeric value that gives calculated p-value from KM survival plot.}
  \item{yTrain}{A Surv object, an output \code{\link[survival:Surv]{Surv}}.}
  \item{survFitOutput}{A survfit object, an output of
                       \code{\link[survival:survfit]{survfit}}.}
  \item{doSkiPlot}{A Boolean value, where \code{TRUE} indicates KM plot is skipped 
                   (not plotted) and vice versa.}
}

\references{
Therneau, T. (2020). A Package for Survival Analysis in R. R package
version 3.1-12, \href{https://CRAN.R-project.org/package=survival}{survival}.

Therneau, T.M., Grambsch, P.M., (2000). Modeling Survival Data: Extending the
Cox Model. \emph{Springer, New York. ISBN 0-387-98784-3}.
}

\seealso{
\code{\link[survival:Surv]{Surv}}, \code{\link[survival:survfit]{survfit}} 
}

\examples{
data(toyCleanedAml)

survival <- toyCleanedAml$survival
inputTime <- setNames(survival[,"Time"]/365, nm=survival[, "PatientID"])
inputEvent <- setNames(survival[,"Dead"], nm=survival[, "PatientID"])
cond <- setNames(survival[,"Risk1"], nm=survival[, "PatientID"])
savePath <- file.path(getwd(), "plotKM.png")

plottedKM <- plotKM(inputTime=inputTime, inputEvent=inputEvent, cond=cond,    
                    plotFile=savePath, verbose=1)
}

\keyword{dplot}
\keyword{survival}
