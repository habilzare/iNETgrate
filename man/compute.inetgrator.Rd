\name{compute.inetgrator}
\alias{compute.inetgrator}

\title{Computes an inetgrator}

\description{It computes the integrator for an arbitrary \code{mu} value. For expert
             users only. Others can use the \code{best.inetgrator} function.}

\usage{
compute.inetgrator(mu, usefuLoci, genePigens, lociPigen, moMe)}

\arguments{
  \item{mu}{
    The numeric value(s) used in network combining in the training dataset. 
    It must be between 0 and 1.}

  \item{usefuLoci}{
    A list of vectors where the names are genes and the values are the corresponding
    loci that contribute to computing eigenloci. This is an output of
    \code{\link{compute.eigenloci}}.}

  \item{genePigens}{
    A list of pigengenes objects with names from
    \code{c("m", "e", "em")} as described in \code{eOrMs} in
    \code{\link{survival.analysis}}.}

  \item{lociPigen}{
    A list of Pigengene object named with genes in the same way as in
    \code{usefuLoci}. This is an output of \code{\link{compute.eigenloci}}.}

  \item{moMe}{
    A character vector of the selected eigengenes to be used 
    e.g., c("ME51.m", "M55.e", "M46.em"). Each eigengene starts with
    "ME" and module number followed by a dot and then the type of data
    used to compute the eigengene. See \code{eOrMs} in
    \code{\link{survival.analysis}}.}
}

\value{
A list with following components:
  \item{mu}{A numeric value used in the network construction of the training dataset. 
    It ranges between \code{0} and \code{1}.}
  \item{pigengenes}{A list of one or more
    \code{\link[Pigengene:pigengene]{pigengene-class}}
    objects, each having the information on genes in the eigengene.
    The names may include "m", "e", and "em" as described
    in \code{eOrMs} in \code{\link{survival.analysis}}.}
  \item{modules}{A list of selected data for each module including the
    weights of genes, the loci for each gene, and if the number of loci is
    more than 1, a Pigengene object for that particular gene.}
}

\seealso{   
\code{\link{compute.eigenloci}}, \code{\link{best.inetgrator}}, 
\code{\link[Pigengene:pigengene]{pigengene-class}} 
}

\examples{
## See the survival.analysis() for an example.
?survival.analysis()
}
\keyword{manip}
