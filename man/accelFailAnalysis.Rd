\name{accelFailAnalysis}
\alias{accelFailAnalysis}

\title{Accelerated Failure Analysis}

\description{
    Performs accelerated failure time model and prediction on each combination 
    of features (e.g., selected modules).}

\usage{
accelFailAnalysis(Data, survival, time2day, eventCol="Dead", riskCol="Risk", 
                  weight=NULL, minRecall4L=0.2, minRecall4H=0.1, until=10, xmin1=0,
                  xmax1=max(survival[, "Time"])*(time2day/365), predType="lp", 
                  doTitle=TRUE, resultPath, 
                  risk2col=c("Low"='green', "Int"='blue', "High"='red'),    
                  abnormalityCol, doAddConfTable=FALSE, favRisk="High", 
                  riskLevel, subSet=NULL, pvalDigits=0, 
                  ylabKm="Survival probability", verbose=0)}

\arguments{
  \item{Data}{
    A matrix of features (eigengenes) values, where rownames are patient IDs and column
    names are features.}

  \item{survival}{
    A matrix or a data frame with minimum of five columns namely: "PatientID", 
    \code{eventCol}, "Time", "Risk1", and "Abnormality". The rownames are the 
    same as "PatientID".}
   
  \item{time2day}{
    A numeric value, that helps calculating survival time in days. E.g. if the "Time" 
    column in survival is in days, time2day=1, else set it to 30 if "Time" is in months.}

  \item{eventCol}{
    A string determining the event column of survival data frame.}

  \item{riskCol}{
    A character vector determining the names of the risk columns of survival data frame.}

  \item{weight}{
    A named numeric vector that determines the weight of each eigengene. Here, 
    names correspond to features and values are weights. If \code{NULL}, all 
    features would have the same weight.}

  \item{minRecall4L}{
    A numeric value of the aimed recall cut--off for low--risk
    predictions, which may not be reached exactly. The higher, the more
    low-risk cases would be identified, but the risk may be relatively higher.}

  \item{minRecall4H}{
    A numeric value of the aimed recall cut--off for high--risk
    predictions, which may not be reached exactly. The higher, the more
    high-risk cases would be identified, but the risk may be relatively lower.}

  \item{until}{
    A numeric value that sets minimum survival time (in years) for an alive case
    to be considered as low risk. Set it to \code{NULL} if you want to create km-plots 
    based on \code{riskCol}.}

  \item{xmax1}{
    A numeric value giving the maximum time (in years) for plotting the K--M curves.}

  \item{xmin1}{
    A numeric value giving the minimum time (in years) for plotting the K--M curves.}

  \item{predType}{
    A character string that takes type input from \code{\link[stats:predict]{predict}}
    function in the stats package.}

  \item{doTitle}{
    Accepts a Boolean value, where \code{TRUE} (default) indicates title will be added
    to plots.}

  \item{resultPath}{
    A character string that gives path to save results.}

  \item{risk2col}{
    A named character vector, that assigns color to the risk groups in KM-plots. 
    Names must be in the \code{riskCol}, and values are colors e.g., 
    risk2col=c(High="red", Low="green", Int="blue").}

  \item{abnormalityCol}{
    A string that determines the column name in survival data that contain risk 
    abnormality information e.g., "Abnormality".}

  \item{doAddConfTable}{
    Accepts a Boolean value, where \code{FALSE} (default) indicates confidence table
    for risk information will not be created.}

  \item{favRisk}{
    A string value that determines a predicted risk group that will be 
    further stratified based on the survival data risk groups. E.g., "High" leads 
    to plotting the survival data risk groups for the subset of data that are
    identified high risk based group on network analysis.}

  \item{riskLevel}{
    A named character vector where values are risk group names (a subset of
    unique(\code{riskCol})), and names are "Low", "Int" and "High". 
    This mapping is necessary for computing the confusion matrix.}

  \item{subSet}{
    A string that determines a risk group from the survival data
    that will be further stratified based on the predicted risk groups 
    e.g., "Int" leads to plotting the predicted risk groups for the subset of 
    data that are determined intermediate risk  group based on survival data.}

  \item{pvalDigits}{
    A numeric value determining the number of decimal digits to be printed on the KM
    plots for p-value.}

  \item{ylabKm}{
    A character vector determining Y--axis label for KM--plot.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
A list with following components:
  \item{subsets}{A list of predicted risk groups information.}
  \item{kmPVals}{A numeric vector of p-values from KM survival plots. }
}
\references{
R Core Team (2019). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
\href{https://www.R-project.org/}{stats}.

Therneau, T. (2020). A Package for Survival Analysis in R. R package
version 3.1-12, \href{https://CRAN.R-project.org/package=survival}{survival}.

Therneau, T.M., Grambsch, P.M., (2000). Modeling Survival Data: Extending the
Cox Model. \emph{Springer, New York. ISBN 0-387-98784-3}.
}

\seealso{
\code{\link[stats:predict]{predict}}, \code{\link{survival.analysis}},
\code{\link{plotKM}}
}

\examples{
data(toyCleanedAml)
data(toyEigengenes)

s1 <- toyCleanedAml$survival
head(s1) ## print

afOut <- accelFailAnalysis(Data=toyEigengenes, survival=s1, 
                           time2day=1, riskCol="Risk1", until=1, xmax1=15, 
                           resultPath=getwd(), abnormalityCol=NULL, 
                           riskLevel=c("Low"="Low", "Int"="Int", "High"="High"), 
                           subSet="Int", verbose=1)
}

\keyword{survival}
