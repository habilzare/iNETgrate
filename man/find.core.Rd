\name{find.core}
\alias{find.core}

\title{Finds the most connected cluster of points}

\description{
    Uses hierarchical clustering and a greedy approach to find "the core" (the most
    connected cluster of loci) for a gene in the given DNA methylation data.
}
\usage{
find.core(Data, Labels, Label1, Label2, verbose=0)}

\arguments{
  \item{Data}{
    A matrix of beta values where column names are probe or loci IDs corresponding 
    to the gene of interest and rownames are patient or case IDs.}

  \item{Labels}{
    A named character vector where names are patient or case IDs and values are
    the lave of the actual risk of patients. See \code{Label1} and  \code{Label2}.}

  \item{Label1}{
    One of the values from \code{Labels} vector that determines a patient group to be used
    to balance the data. E.g. low--risk group defined as "Low" in
    toyCleanedAml.}

 \item{Label2}{
    One of the values from \code{Labels} vector that determines a patient group to be used
    to balance the data. E.g. high--risk group defined as "High" in
    toyCleanedAml.}

  \item{verbose}{
    The integer level of verbosity. \code{0} means silent, and higher values produce 
    more details of the computation.}

}

\value{
A list with following components:
  \item{PC1}{The matrix of inferred (projected) eigengenes where rows are
             samples and columns are modules.}
  \item{pigenObject}{An object of Pigenengene class.}
  \item{lociNames}{A vector of loci that we used for computing methylation level 
                   of the related gene (using PCA, a weighted average value).}
  \item{aveWeight}{Average similarity level of the loci group that relates to a gene}
  \item{communities}{Clustered loci}
  \item{num}{Number of unique loci clusters}
}
\references{
Amir Foroushani \emph{et al.} (2016) Large-scale gene network analysis reveals the significance
of extracellular matrix pathway and homeobox genes in acute myeloid leukemia: an
introduction to the Pigengene package and its applications, \emph{BMC MedicalGenomics}.
}

\examples{
data(toyCleanedAml)

unionGenes <- c("BAI1", "CAPNS1", "CHD1", "CYC1", "EXT1", "FYB", "GFRA1", "HLX", 
                "IDH3B", "PSMD8", "SGTA", "SLC1A5", "SSFA2", "TFAM", "TGIF1", 
                "UBC", "ZNF70", "ST8SIA4", "ZNF485", "FOXN3", "HCG4", "FOXK1", 
                "C2orf44")

l2g <- toyCleanedAml$locus2gene
survival <- toyCleanedAml$survival

loci <- l2g[is.element(l2g[, "Gene_Symbol"], unionGenes), "probeID"]
patientLabel <- setNames(as.character(survival$Risk1), rownames(survival))
inBoth <- intersect(colnames(toyCleanedAml$dnam), names(patientLabel))
LabelsIn <- patientLabel[is.element(names(patientLabel), inBoth)]
dnamData <- t(toyCleanedAml$dnam[loci, inBoth])

core <- find.core(Data=dnamData, Labels=LabelsIn, Label1="Low", Label2="High", 
                  verbose=1)
}

\keyword{graphs}
\keyword{cluster}
