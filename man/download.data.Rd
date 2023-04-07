\name{download.data}
\alias{download.data}

\title{Download and save TCGA data}
\description{
    Uses TCGAbiolinks to query, download and save TCGA open-access data as an R
    object. Legacy data is downloaded whenever available.}

\usage{
download.data(dataProject, savePath, doDnam=TRUE,  doExpr=TRUE, 
              doClinical=TRUE, doMirna=FALSE, legacy=TRUE, verbose=0)}

\arguments{
  \item{dataProject}{
    A character vector indicating the name for project as listed in TCGA portal.}

  \item{savePath}{
    A character vector that provides path to the folder where data is saved.}

  \item{doDnam}{
    Accepts a Boolean value, where \code{TRUE} (default) indicates DNA methylation data 
    will be downloaded.}

  \item{doClinical}{
    Accepts Boolean value, where \code{TRUE} (default) indicates clinical data will be downloaded.}

  \item{doExpr}{
    Accepts Boolean value, where \code{TRUE} (default) indicates gene expression data
    will be downloaded.}

  \item{doMirna}{
    Accepts a Boolean value, where \code{TRUE} indicates micro RNA data will be downloaded.
    Default is set to \code{FALSE}.}

  \item{legacy}{
    Accepts a Boolean value, where \code{TRUE} (default) indicates legacy data will be 
    downloaded.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}
}

\value{
A list with data saved as R object, paths to saved data, and parameters used.
Some of the list elements are:
  \item{clinical}{A data frame of clinical data where rows are patient IDs;
                  columns are clinical parameters.}

  \item{genExpr}{A data frame of gene expression data where rows are genes and 
                 columns are sample IDs.}

  \item{rawDnam}{A SummarizedExperiment object with DNA methylation data.}

  \item{genExprSampleInfo}{A data frame where rows are samples and
                           columns are associated sample information.}
}

\references{
Colaprico, A., Silva T.C., Olsen, C., Garofano, L., Cava, C., Garolini, D.,
Sabedot, T.,  Malta, T., Pagnotta, S.M., Castiglioni, I., Ceccarelli, M.,
Houtan Noushmehr, G.B.,  (05 May 2016) TCGAbiolinks: An R/Bioconductor package 
for integrative analysis of TCGA data \emph{Nucleic Acids Research} \bold{44(8): e71}.
(doi:10.1093/nar/gkv1507); 
\href{https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html}{TCGAbiolinks}

Silva, T. C. \emph{et al.} (2016). TCGA Workflow: Analyze cancer genomics and 
epigenomics data using Bioconductor packages, \emph{F1000Research} \bold{5}.

Mounir, M. \emph{et al.} (2019). New functionalities in the TCGAbiolinks package 
for the study and integration of cancer data from GDC and GTEx, 
\emph{PLoS computational biology} \bold{15.3 : e1006701}.
}

\seealso{
 \code{\link[TCGAbiolinks:GDCquery]{GDCquery}},
 \code{\link[TCGAbiolinks:GDCdownload]{GDCdownload}},
 \code{\link[TCGAbiolinks:GDCprepare]{GDCprepare}}, \code{\link{clean.allData}} 
}

\examples{
\dontrun{ ## New TCGAbiolinks is buggy (2.24.3 < Version  <= 2.27.2).  
tcgaAml <- download.data(dataProject="TCGA-LAML", savePath=getwd(), 
                         doDnam=FALSE, doClinical=TRUE, 
                         doExpr=FALSE, doMirna = FALSE, legacy=TRUE)
}
}

\keyword{data}
