\name{inferEigengenes}
\alias{inferEigengenes}

\title{Infers features for the test(validation) dataset}

\description{
    It computes features (eigengenes) for the test (validation) dataset.}

\usage{
inferEigengenes(inetgrator, dnam=NULL, expr, verbose=0)}

\arguments{
  \item{inetgrator}{
    A list with at least 3 components namely "mu", "pigengenes", and
    "modules" similar to the output from \code{\link{bestInetgrator}}.}

  \item{dnam}{
    A matrix of beta values where rows are samples and columns are loci.}

  \item{expr}{
    A matrix of expression values where rows are samples and columns are genes.}

  \item{verbose}{
    An integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
A list with following components:
  \item{inetgrator}{A list same as the input \code{inetgrator}.}
  \item{eigengeneS}{A matrix with samples on rows, the inferred features
                    (eigengenes) on columns.}
  \item{missingInExpr}{A character vector of genes missing in the test dataset.}
}

\seealso{
\code{\link{bestInetgrator}}, \code{\link{computeEigenloci}}
}

\examples{
## See the analyzeSurvival() for an example.
?analyzeSurvival()
}

\keyword{manip}
