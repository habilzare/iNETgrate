\name{iNETgrate}
\alias{iNETgrate}

\title{iNETgrate - Runs entire iNETgrate pipeline.}

\description{
    Runs the entire iNETgrate pipeline; from cleaning the data to finding the best
    modules and performing survival analysis.
    Specifically, it computes the eigenloci
    (weighted average of beta values per gene), identifies the gene modules 
    using coexpression network analysis, computes eigengenes, identifies best
    modules, reclassifies patients into risk groups based on the
    eigengene values, and performs survival analysis.}

\usage{
iNETgrate(Data, clinSettings, mus=(0:10)/10, saveDir="iNETgrate", 
          annLib="Auto", isFfpe=c(FALSE), doRemoveTOM=TRUE, 
          minModuleSize=5, corMethod="pearson", doReturNetworks=FALSE, 
          RsquaredCut=0.75, combiningMu=NA, favRisk="High", subSet="Int", 
          xmax1=15, xmin1=0, doCox=TRUE, eOrMs="e", time2day=1, 
          until=1, minRecall4L=0.2, minRecall4H=0.05, verbose=0)
}

\arguments{

  \item{Data}{
    A list with minimum of four elements named \code{rawDnam}, 
    \code{clinical}, \code{genExpr}, and \code{genExprSampleInfo}.
    Here, \code{rawDnam} is an object
    of class \code{SummarizedExperiment}, \code{genExpr} is the matrix
    of gene expression values with genes on rows and patients on columns, and \code{clinical}
    and \code{genExprSampleInfo} are data frames with patients on rows.
    Example: \code{\link{toyRawAml}}.}

  \item{clinSettings}{
    A vector where names are "patientIDCol", "eventCol", "timeCol", "riskCatCol", 
    "riskFactorCol", "otherLabel", "event", "riskHigh", and "riskLow", and the 
    values are the corresponding information from the clinical data. In particular, 
    values for a name ending in "Col" must be the corresponding column name in
    the clinical dataframe e.g., "patientIDCol" in the name of the
    column that has the patient ID. "event" is the event of interest
    (e.g., "Dead" may mean dead due to the disease) as indicated in
    "eventCol". "riskHigh", and "riskLow" are the values of "riskCatCol"
    corresponding to high-- and low--risk, respectively. 
    See the \code{\link{clean.allData}} function and  the example for
    further details.}

  \item{mus}{
    A numeric vector assigning mu value(s) in the range [0,1].
    A higher value for mu will  result in more emphasize on the DNA
    methylation data than gene expression data. In particular, mu=0
    would lead to ignoring DNA methylation data, with mu=1 gene
    expression data would have no impact on the network, and with mu=0.5
    DNA methylation and gene expression would have equal contribution to
    the network.}

  \item{saveDir}{
    A character string describing the path to directory for saving plots and output.}

  \item{annLib}{
    A character string that determines the annotation (reference) library to be
    used. The default is "IlluminaHumanMethylation450kanno.ilmn12.hg19".}

  \item{isFfpe}{
    A character vector that takes Boolean values to indicate inclusion or
    exclusion of FFPE sample. c(FALSE, TRUE) will include all samples, TRUE 
    will include only FFPE samples, and FALSE will exclude FFPE samples.
    Set it to \code{NULL} to disable it.}

  \item{doRemoveTOM}{
    A Boolean determining if TOM files must removed. Default is set to
    \code{TRUE}.}
  
  \item{minModuleSize}{
    A numeric value that specifies the minimum number of genes per module.}
  
  \item{corMethod}{
    A string that specifies the correlation method to be used. Accepted
    options include "spearman" or "pearson", and default is set to "pearson".}

  \item{doReturNetworks}{
    A Boolean value determining whether to return \code{exprNetwork} and \code{dnamNetwork}, 
    which are relatively big objects (typically Gbs). Default: \code{FALSE}.}

  \item{RsquaredCut}{
    A threshold in the range [0,1] used to estimate the power (beta
    value) for soft thresholding by WGCNA. A higher value can increase
    the power to which the similarity matrix is raised. For technical
    use only. 
    See  \code{\link[WGCNA:pickSoftThreshold]{pickSoftThreshold}} for more details.}

  \item{combiningMu}{
    A numeric vector determining the mu value(s) in the range [0,1]
    used for eigengene computation.  If set to NA, the same value that
    was used for network construction will be used for eigengene
    computation. See \code{mus}.}

  \item{favRisk}{
    The value must be one of these three: "High", "Int" or "Low".
    Determines the risk category predicted by network classifier that shall be compared to the 
    input \code{Labels} risk categories.}

  \item{subSet}{
    The value must be one of these three: "High", "Int" or "Low".
    Determines the risk category in \code{Labels} that shall be compared to
    the risk category predicted by network classifier.}

  \item{xmax1}{
    A numeric value determining the maximum time (in years) for plotting the K--M curves.}

  \item{xmin1}{
    A numeric value determining the minimum time (in years) for plotting the K--M curves.}

  \item{doCox}{
    Accepts a Boolean value, where \code{TRUE} (default) indicates Cox analysis will be
    performed for module selection.}

  \item{eOrMs}{
    A character vector to determine which eigengenes (stored at netPath) 
    to be used as input data for survival analysis, e.g., c("m", "e",
    "em"), where "e" means "only gene expression", "m" means
    "DNA methylation", and "em" means the combination of both data types
    is used to compute an eigengene for a module.}

  \item{time2day}{
    A numeric value used to calculate survival time (in days). E.g. if the "Time" 
    column in \code{survival} data frame is in days, time2day=1, set it to 
    30 if "Time" is in months, etc.}

  \item{until}{
    A numeric value describing minimum survival time (in years) for a
    case to be considered as low risk, e.g., if an AML patient was alive for at 
    least 2 years after diagnosis AND was alive at the last followup
    time, then we can consider this case low--risk.}

  \item{minRecall4L}{
    A numeric value of the aimed recall cut--off for low--risk
    predictions, which may not be reached exactly. The higher, the more
    low--risk cases would be identified, but the risk in the group may
    be relatively higher. Do not change the default if you are not an
    expert.}

  \item{minRecall4H}{
    A numeric value of the aimed recall cut--off for high--risk
    predictions, which may not be reached exactly. The higher, the more
    high-risk cases would be identified, but the risk in the group may
    be relatively lower. Do not change the default if you are not an
    expert.}

  \item{verbose}{
    An integer level of verbosity: \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\details{
    This is a wrapper function to run the iNETgrate pipeline in several
    steps. First, the input data is cleaned by dropping any missing data and
    imputing data as needed. Next, we exclude the genes for which there is no or negligible
    correlation between expression levels nor DNA methylation and
    survival data. We compute a weighted average of  beta values per
    gene, and build a network using both gene expression as well as DNA methylation
    data. Gene modules are identified based on coexpression and DNA methylation 
    correlation between genes. Then, the eigengenes for each module 
    and each sample are calculated, where the expression of an eigengene of a
    module in a sample is the weighted average of the expression of the genes in 
    that module in that sample. Technically, an eigengene is the first principal 
    component of the gene expression levels in a module. Using PCA, we
    ensures that the maximum variance across all the training samples
    is explained by the eigengene. Next, we select a combination of 3
    modules or less, based on a Cox regression analysis of survival data
    and the corresponding module eigengenes. Along the way, several
    self explanatory directories, and plots are created and stored under
    \code{saveDir}.}


\value{
A list with following components:
  \item{call}{The call that created the results.}
  \item{unionGenes}{A vector of gene names that are selected to make the
    network.}
  \item{Labels}{A vector where names are patient IDs and values are
                the corresponding risk levels.}
  \item{eigenloci}{A matrix of DNA methylation level per genes (columns)
    with each row corresponding to a patient.}
  \item{bestInet}{A list of outputs from \code{\link{best.inetgrator}} function.}
}

\references{
R Core Team (2019). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria.
\href{https://www.R-project.org/}{stats}.

Foroushani, Amir, et al. (2017) Large-scale gene network analysis
reveals the significance of extracellular matrix pathway and homeobox
genes in acute myeloid leukemia: an introduction to the Pigengene
package and its applications, \emph{BMC medical genomics}, \bold{10.1},
doi.org/10.1186/s12920-017-0253-6.

Langfelder, P. and Horvath, S., (2008) WGCNA: an R package for weighted correlation 
network analysis, \emph{BMC Bioinformatics}, \bold{9:559},
doi:10.1186/1471-2105-9-559.

Langfelder P., Horvath S. (2012) Fast R Functions for Robust Correlations and
Hierarchical Clustering. \emph{Journal of Statistical Software}, \bold{46(11),1-17}. 
\href{http://www.jstatsoft.org/v46/i11/}{WGCNA}.
}

\seealso{
\code{\link[WGCNA:pickSoftThreshold]{pickSoftThreshold}},
\code{\link{clean.allData}}, \code{\link{elect.genes}},
\code{\link{make.network}}, \code{\link{combined.eigengenes}},
\code{\link{best.inetgrator}}, \code{\link{survival.analysis}},
\code{\link{toyRawAml}}
}

\note{
The individual functions are exported to facilitate running the pipeline 
step--by--step and in a customized way.}

\examples{

## Preparing data:
data(toyRawAml)

## Set up clinical settings
clinSettings <- c("patientIDCol"="bcr_patient_barcode",
"eventCol"="vital_status",
"timeCol"="days_to_last_followup",
"riskCatCol"="acute_myeloid_leukemia_calgb_cytogenetics_risk_category",
"riskFactorCol"="cytogenetic_abnormalities",
"otherLabel"="NULL",
"event"="Dead",
"riskHigh"="Poor",
"riskLow"="Favorable")
print(clinSettings)

## Printing few rows and columns of clinical data
head(toyRawAml$clinical[ , clinSettings[c("patientIDCol","eventCol", "timeCol", 
                                          "riskCatCol", "riskFactorCol")]])

inetgrateRes <- iNETgrate(Data=toyRawAml, clinSettings=clinSettings, mus=0.6)
}

\keyword{optimize}
\keyword{cluster}
\keyword{classif}
\keyword{models}
\keyword{survival} 
