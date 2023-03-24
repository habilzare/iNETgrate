\name{iNETgrate-package}

\alias{iNETgrate-package}

\docType{package}
\title{Integrates DNA methylation data with gene expression in a single gene network}

\description{
The iNETgrate package can be used to build a correlation network in which nodes
are genes. DNA methylation and gene expression data are integrated to define the
connections between genes. This network is used to identify modules (clusters) of
genes. The biological information in each of the resulting modules is represented by an
eigengene. These biological signatures can be used as features e.g., for classification
of patients into risk categories. The resulting biological signatures are very robust
and give a holistic view of the underlying molecular changes.}  

\details{
\tabular{ll}{
Package: \tab iNETgrate\cr
Type: \tab Package\cr
Version: \tab 0.5.24\cr
Date: \tab 2021-08-25\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Isha Mehta, Ghazal Ebrahimi, Hanie Samimi, and Habil Zare

Maintainer: Habil Zare <zare@uthscsa.edu>
}

%% \references{}

\seealso{
\code{\link[iNETgrate]{iNETgrate}}, \code{\link[Pigengene]{Pigengene-package}}, 
\code{\link[TCGAbiolinks:GDCdownload]{GDCdownload}},
\code{WGCNA::\link[WGCNA]{blockwiseModules}}
}

\examples{
%% See iNETgrate().
?iNETgrate
}

\keyword{package}
\keyword{datasets}
\keyword{documentation}
\keyword{models}
