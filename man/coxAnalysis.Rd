\name{coxAnalysis}
\alias{coxAnalysis}

\title{Cox analysis}

\description{Function to perform Cox regression analysis.}

\usage{
coxAnalysis(dataCategory="training", inputTime, inputEvent, inputMatrix, 
            outPath, weight=NULL, verbose=0)
}

\arguments{
  \item{dataCategory}{
    A string describing the type of data analyzed e.g., "Training
    data", "validation data", or "all data".}

  \item{inputTime}{
    A named numeric vector where names are patient IDs and values are survival 
    time (in years).}

  \item{inputEvent}{
    A named binary vector where names are patient IDs and values are \code{0}
    and \code{1}. Here, \code{1} indicates the occurance of event and vice versa.}

  \item{inputMatrix}{
    A matrix of projected eigengenes for each sample, where row names are case
    IDs and column names are modules.}

  \item{outPath}{
    Output path for coxAnalysis.}

  \item{weight}{
    Same as \code{weights} argument in \code{\link[glmnet:glmnet]{glmnet}}.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}

}
\value{
A character vector of ranked input features (eigengenes) from high to low based 
on Cox analysis.}

\references{
Friedman, J., Hastie, T., Tibshirani, R. (2010). Regularization Paths for Generalized 
Linear Models via Coordinate Descent; \emph{Journal of Statistical Software}
\bold{33(1), 1-22}; \href{http://www.jstatsoft.org/v33/i01/}{glmnet}.

Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011). Regularization Paths 
for Cox's Proportional Hazards Model via Coordinate Descent; \emph{Journal of Statistical 
Software} \bold{39(5), 1-13}; \href{http://www.jstatsoft.org/v39/i05/}{glmnet}.

Therneau, T. (2020). A Package for Survival Analysis in R; \emph{R package
version 3.1-12}; \href{https://CRAN.R-project.org/package=survival}{survival}.

Therneau, T.M., Grambsch, P.M., (2000). Modeling Survival Data: Extending the
Cox Model; \emph{Springer, New York} \bold{ISBN 0-387-98784-3}.
}

\seealso{
\code{\link[glmnet:glmnet]{glmnet}}
}

\examples{
data(toyCleanedAml)
data(toyEigengenes)
survival <- toyCleanedAml$survival
inputTime <- setNames(survival[,"Time"]/365, nm=survival[,"PatientID"])
inputTime <- inputTime[inputTime > 0]
inBoth <- intersect(names(inputTime), rownames(toyEigengenes))
inputTime <- inputTime[is.element(names(inputTime), inBoth)]
inputEvent <- setNames(survival[is.element(survival$PatientID, inBoth), "Dead"], 
                       nm=names(inputTime))

coxOut <- coxAnalysis(dataCategory="allSamples", inputTime=inputTime, 
                      inputEvent=inputEvent, inputMatrix=toyEigengenes[inBoth,],
                      outPath=tempdir(), verbose=1)
}

\keyword{glmnet}
\keyword{Surv}
