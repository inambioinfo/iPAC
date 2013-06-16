\name{iPAC-package}
\alias{iPAC-package}
\alias{iPAC}
\docType{package}
\title{
Identification of Protein Amino acid mutation Clustering
}
\description{
This package finds mutation clusters on the amino acid level while taking into acount the protein structure.
}
\details{
The \emph{iPAC} package finds clusters of non-synonomous amino acid mutations via a two step process. First, \emph{iPAC} links the amino acid positions in 3D space with their ``canonical" ordering. The canonical ordering is considered to be the sequence of amino acids when the protein is viewed in linear form. Next, the canonical protein ordering is used to match the 3D positional information with the mutation file. The mutation file is a matrix of 0's (no mutation) and 1's (mutation) where each column represents an amino acid and each row represents an individual sample (test subject, cell line, etc). Thus if column i in row j had a 1, that would mean that the ith amino acid for person j had a nonsynonomous mutation. The last pre-processing step then maps the amino acids onto a 1D space while preserving (as best as possible) local neighbor relationships. The NMC algorithm (see \emph{Ye. et. al.}) is then executed to identify the significant mutational clusters. Finally, the information is mapped back into the original 3D space and reported back to the user.
}
\author{
Author: Gregory A. Ryslik, Hongyu Zhao

Maintainer: Gregory A. Ryslik <gregory.ryslik@yale.edu>

}
\references{
%We would like to extend a special thank you to Dr. Jingjing Ye and her group for their assistance with regards to the initial NMC algorithm.

Ye et. al., Statistical method on nonrandom clustering with application to somatic mutations in cancer. \emph{BMC Bioinformatics}. 2010. doi:10.1186/1471-2105-11-11.

Bioconductor: Open software development for computational biology and bioinformatics R. Gentleman, V. J. Carey, D. M. Bates, B.Bolstad, M. Dettling, S.
  Dudoit, B. Ellis, L. Gautier, Y. Ge, and others 2004, Genome Biology, Vol. 5, R80
  
Borg, I., Groenen, P. (2005). Modern Multidimensional Scaling: theory and applications (2nd ed.). New York: Springer-Verlag.
}

\examples{
\dontrun{
#Extract the data from a CIF file and match it up with the canonical protein sequence.
#Here we use the 3GFT structure from the PDB, which corresponds to the KRAS protein.
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions<-get.AlignedPositions(CIF,Fasta, "A")

#Load the mutational data for KRAS. Here the mutational data was obtained from the
#COSMIC database (version 58). 
data(KRAS.Mutations)

#Identify and report the clusters. 
ClusterFind(mutation.data=KRAS.Mutations, 
							position.data=KRAS.Positions$Positions,
							create.map = "Y",Show.Graph = "Y")
}

}

\keyword{package}
