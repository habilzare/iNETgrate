\name{computeEigengenes}
\alias{computeEigengenes}

\title{Compute eigengenes of combined modules.}
\description{
    Computes eigengenes (features) for network modules obtained from expression and
    dnam data combined.}

\usage{
computeEigengenes(genExpr=NULL, eigenloci=NULL, netPath, geNames,
                    Labels, Label1, Label2, mus, combiningMu=NA,
                    survival, event="Dead", doIgnoreNas=FALSE, verbose=0, 
                    mu2modules)
}

\arguments{
  \item{genExpr}{
    A matrix of gene expression values where row names are gene IDs, and column 
    names are patient IDs.}

  \item{eigenloci}{
    A matrix of eigenloci where row names are patient IDs, and column names are 
    gene IDs.}

  \item{netPath}{
    A string that determines the path to a folder where output object
    \code{combinedNetwork} is saved.}

  \item{geNames}{
    A character vector of selected genes that become the nodes of the network. Both
    gene expression and eigenloci data must be available for these genes e.g., 
    unionGenes (output of \code{\link{electGenes}}).}

  \item{Labels}{
    A named character vector that determines the risk group for each patient. Here,
    the names are the patient IDs and the values are the risk group labels. It is needed for
    sample balancing when computing the eigenloci.}

  \item{Label1}{
    One of the values from \code{Labels} vector that defines a patient group to be used to
    balance data e.g., low--risk group defined as "Favorable" in toyAmlData.}

 \item{Label2}{
   One of the values from \code{Labels} vector that defines a patient group to be used to 
   balance data e.g., high--risk group defined as "Poor" in toyAmlData.}

  \item{mus}{
    A numeric value or vector. Lambda values that are used for network construction.}

  \item{combiningMu}{
    A numeric vector determining the mu value(s) used for eigengene computation. 
    If set to NA, the same value that was used for network construction, will be 
    used for eigengene computation.}

  \item{survival}{
    A data frame where row names are patient IDs and must have minimum of four columns
    namely: "Time" (to event), \code{event}, "Abnormality" and "Risk1".}

  \item{event}{
    A string determining the name of the \code{event} column in the survival data. 
    It stores binary values 1 and 0 that correspond to the occurance or non-occurance 
    of the event respectively.}

  \item{doIgnoreNas}{
    A Boolean with a default value of \code{False}.}

  \item{verbose}{
    An integer level of verbosity. 0 means silent, and higher values produce more 
    details of the computation.}

  \item{mu2modules}{
    A matrix where rows are mus and values are module information. It is an 
    output of \code{\link{makeNetwork}} function.}
}

\value{
A list with following components:
  \item{patients}{A character vector of patientIDs that are included in
                   eigengene calculations. }
  \item{lcombine}{A numeric vector where names are "e", "m", and "em". The
                   value corresponding to "em" is same as mus value. If
                   genExpr and eigenLoci inputs are available, "e" and "m" are
                   set to 0 and 1 respectively.}
  \item{eigengeneFiles}{A matrix where rows are mus, columns are
                        lcombine, and values are paths to eigengene RData files.}
}
\references{
Foroushani, A. \emph{et al.} (2016) Large-scale gene network analysis reveals the significance
of extracellular matrix pathway and homeobox genes in acute myeloid leukemia: an
introduction to the Pigengene package and its applications, \emph{BMC MedicalGenomics}.
}

\seealso{
\code{\link{makeNetwork}}
\code{\link[Pigengene:compute.pigengene]{compute.pigengene}}, 
\code{\link[Pigengene:project.eigen]{project.eigen}}
}

\examples{
data(toyCleanedAml)
data(toyComputEloci)
eigenloci <- toyComputEloci$eigenloci 
patientLabel <- as.matrix(toyCleanedAml$survival)[,"Risk1"]

# Path to the network
netPath <- file.path(tempdir(), "net")
dir.create(netPath, showWarnings=FALSE)
madeNetwork <- makeNetwork(genExpr=toyCleanedAml$genExpr, eigenloci=eigenloci,
                            geNames=colnames(eigenloci), mus=c(0.6), 
                            outPath=netPath, verbose=0)

eigengenes <- computeEigengenes(genExpr=toyCleanedAml$genExpr, 
                                  eigenloci=eigenloci, netPath=netPath, 
                                  geNames=colnames(eigenloci), Labels=patientLabel, 
                                  Label1="Low", Label2="High", 
                                  mus=c(0.6), survival=toyCleanedAml$survival,
                                  mu2modules=madeNetwork$mu2modules)
}
\keyword{multivariate}
