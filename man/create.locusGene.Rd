\name{create.locusGene}
\alias{create.locusGene}

\title{Mapping one locus to one gene.}

\description{
Creates a dataframe that maps each loci to its corresponding gene and vice versa.}

\usage{
create.locusGene(locus2gene="Auto", genesColName="Gene_Symbol", lociColName="probeID", 
                 verbose=0)}

\arguments{
  \item{locus2gene}{
    A dataframe with minimum 4 columns including \code{genesColName},
    \code{lociColName}, "Chromosome", and "Genomic_Coodinate". Default "Auto" 
    gets loci annotations from "IlluminaHumanMethylation450kanno.ilmn12.hg19".}

  \item{genesColName}{
    A string indicating the column name in the \code{locus2gene} dataframe that
    contain genes information.}

  \item{lociColName}{
    A string indicating the column name in the \code{locus2gene} dataframe that
    contain loci information.}

  \item{verbose}{
    The integer level of verbosity. 0 means silent, and higher values produce 
    more details of the computation.}
}

\value{
  A list with following components:
    \item{noGenes}{A character vector of loci names that lack gene assignment.}
    \item{locus2oneGene}{A dataframe with one locus mapped to one gene.}
}

\examples{
data(toyCleanedAml)

locus2oneGene <- create.locusGene(locus2gene=toyCleanedAml$locus2gene, 
                                  verbose=1)
}
\keyword{manip}
\keyword{datasets}
