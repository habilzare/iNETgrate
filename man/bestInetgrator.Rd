\name{bestInetgrator}
\alias{bestInetgrator}

\title{Identifies the bestInetgrator}

\description{
    It creates the best integrator object by collecting all the weights
    corresponding to the best \code{mu} value from the results of the
    training dataset. The inetgrator can then be used to infer the
    eigengenes in validation datasets using the
    \code{\link{inferEigengenes}} function.}

\usage{
bestInetgrator(bestPvalues, usefuLoci, lociPigen, netPath, verbose=0)}

\arguments{
  \item{bestPvalues}{
    A data frame with at least three columns named: "pvalue", "mu", and "features". 
    This is similar to the \code{bestPvalues} output from \code{\link{analyzeSurvival}}.}

  \item{usefuLoci}{
    A list of vectors where the names are genes and the values are the
    corresponding loci that contribute to computing eigenloci. This is
    an output of \code{\link{computEigenloci}}.}

  \item{lociPigen}{
    A list of genes (same as \code{usefuLoci}'s genes and their related 
    Pigengene object (based on their mapped loci). This is an output of
    \code{\link{computEigenloci}}.}

  \item{netPath}{
    A string determining the path to the network folder, where the pigengene and 
    eigengene objects saved for selected modules of each mu.}

  \item{verbose}{
    An integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
   The output of the \code{\link{computeInetgrator}} function
   corresponding to the best \code{mu} based on the analyzeSurvival.
}

\seealso{
  \code{\link{inferEigengenes}},
  \code{\link{analyzeSurvival}}, \code{\link{computEigenloci}},
  \code{\link{computeInetgrator}},
  \code{\link[Pigengene:pigengene]{pigengene-class}}
}

\examples{
## See the analyzeSurvival() for an example.
?analyzeSurvival()
}

\keyword{manip}
