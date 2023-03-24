\name{make.network}
\alias{make.network}

\title{Make network}

\description{
  Uses the adjacency to create a network from input data. Then, uses the
  mu value to create a combined network from the expression network and
  DNAm network.}

\usage{
make.network(genExpr, eigenloci, geNames, mus, doRemoveTOM=TRUE, outPath, 
             minModuleSize=5, corMethod="pearson", doReturNetworks=FALSE, 
             RsquaredCut=0.75, doSaveCombined=FALSE, verbose=0)}

\arguments{
  \item{genExpr}{
    A gene expression matrix where row names are gene IDs and column names are 
    patient IDs.}

  \item{eigenloci}{
    An eigenloci matrix where row names are patient IDs and column names are 
    gene IDs.}

  \item{geNames}{
    A character vector of selected genes that become the nodes of the network. 
    Both gene expression and eigenloci data must be available for these genes
    (e.g., \code{unionGenes} in the output of \code{\link{elect.genes}}).}

  \item{mus}{
    A numeric vector assigning mu value(s).}

  \item{doRemoveTOM}{
    A Boolean determining if TOM files must removed. Default is set to
    \code{TRUE}. For technical users only.}
  
  \item{outPath}{
    A string that determines the path to the folder where plots and output 
    will be saved.}

  \item{minModuleSize}{
    A numeric value that specifies the minimum number of genes per
    module. See \code{\link[WGCNA:blockwiseModules]{blockwiseModules}}.}
  
  \item{corMethod}{
    A string that specifies the correlation method to be used. Accepted
    options include "spearman" or "pearson", and default set to "pearson".}

  \item{doReturNetworks}{
    A Boolean value determining whether to return \code{exprNetwork} and \code{dnamNetwork}, 
    which are relatively big objects (typically Gbs). Default \code{FALSE}.}

  \item{RsquaredCut}{
    A threshold in the range [0,1] used to estimate the power. A higher value can 
    increase power for technical use only. See 
    \code{\link[WGCNA:pickSoftThreshold]{pickSoftThreshold}} for more details.}

  \item{doSaveCombined}{
    If \code{TRUE}, the results from \code{\link[Pigengene:combine.networks]{combine.networks}}
    function will be saved on storage. Leave it as \code{FALSE} not to waste typically 1-2
    GBs disk space.}

  \item{verbose}{
    An integer level of verbosity. \code{0} means silent, and higher values produce more 
    details of the computation.}
}

\value{
A list with following components:
  \item{corMethod}{A character value same as the input.}
  \item{mus}{A numeric vector of mu value(s) same as the input.}
  \item{mu2modules}{A matrix where rows are mus and values are module
                    information that is an output of 
                    \code{\link[Pigengene:combine.networks]{combine.networks}}.}
  \item{powerVector}{A numeric vector of power values named with the corresponding mus. 
                     See \code{\link[WGCNA:pickSoftThreshold]{pickSoftThreshold}} documentation.}
  \item{outliersNumber}{A named numeric vector where names are mu value(s) and
                        values are number of outlier genes (genes in module 0).}

If doReturNetworks is TRUE the following components are returned:
  \item{exprNetwork}{An adjacency matrix of gene expression.}
  \item{dnamNetwork}{An adjacency matrix of DNA methylation.}
}
  
\references{
Langfelder, P. and Horvath, S., (2008) WGCNA: an R package for weighted correlation 
network analysis, \emph{BMC Bioinformatics}, \bold{9:559}, doi:10.1186/1471-2105-9-559

Langfelder P., Horvath S. (2012) Fast R Functions for Robust Correlations and
Hierarchical Clustering. \emph{Journal of Statistical Software}, \bold{46(11),1-17}. 
\href{http://www.jstatsoft.org/v46/i11/}{WGCNA}

Foroushani, A. \emph{et al.} (2016) Large-scale gene network analysis reveals the significance
of extracellular matrix pathway and homeobox genes in acute myeloid leukemia: an
introduction to the Pigengene package and its applications, \emph{BMC MedicalGenomics}.
}

\seealso{
\code{\link[Pigengene:combine.networks]{combine.networks}},
\code{\link[WGCNA:pickSoftThreshold]{pickSoftThreshold}},
\code{\link[WGCNA:blockwiseModules]{blockwiseModules}}, \code{\link{elect.genes}}
}

\examples{
data(toyCleanedAml)
data(toyComputEloci)
eigenloci <- toyComputEloci$eigenloci

netPath <- file.path(getwd(), "net")
dir.create(netPath, showWarnings=FALSE)
madeNetwork <- make.network(genExpr=toyCleanedAml$genExpr, eigenloci=eigenloci,
                            geNames=colnames(eigenloci), mus=c(0.6), 
                            doRemoveTOM=FALSE, outPath=netPath, verbose=1)
}

\keyword{graph}
\keyword{cluster}
