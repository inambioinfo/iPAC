\name{get.Remapped.Order}
\alias{get.Remapped.Order}
\title{
get.Remapped.Order
}
\description{
get.Remapped.Order returns the reordering of the culled mutation matrix using the remapper of choice (linear vs MDS).
}
\usage{
get.Remapped.Order(mutation.data, position.data, method = "MDS", OriginX = min(position.data[,4]), OriginY=min(position.data[,5]), OriginZ=min(position.data[,6]))
}
\arguments{
  \item{mutation.data}{
	A matrix of 0's (no mutation) and 1's (mutation) where each column represents an amino acid in the protein and each row represents an individual sample (test subject, cell line, etc). Thus if column i in row j had a 1, that would mean that the ith amino acid for person j had a nonsynonomous mutation.
}
  \item{position.data}{
	A dataframe consisting of five columns: 1) Residue Name, 2) Amino Acid number in the protein, 3) Side Chain, 4) X-coordinate, 5) Y-coordinate and 6) Z-coordinate. Please see \emph{get.Positions} and \emph{get.AlignedPositions} for further information on how to construct this matrix. 
}
  \item{method}{
	You can select whether you want a "MDS" or "Linear" approach in order to map the protein into a 1D space.
}
  \item{OriginX}{
	If the "Linear" method is chosen, this specifies the x-coordinate part of the fixed point.
}
  \item{OriginY}{
	If the "Linear" method is chosen, this specifies the y-coordinate part of the fixed point.
}
  \item{OriginZ}{
	If the "Linear" method is chosen, this specifies the z-coordinate part of the fixed point.
}
}
\details{
  The culled mutation matrix is the mutation matrix after the amino acids that have missing positional data have been removed. The amino acid positions that have missing positional data will not be displayed in the result.
 
}

\note{
The returned value shows the remapped order of the culled mutation matrix.

This method is still in beta. If any bugs found, please email the package authors.
}

\examples{
#Extract the data from a CIF file and match it up with the canonical protein sequence.
#Here we use the 3GFT structure from the PDB, which corresponds to the KRAS protein.
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.sanger.ac.uk/perl/genetics/CGP/cosmic_sequence?ln=KRAS;format=aa"
KRAS.Positions<-get.AlignedPositions(CIF,Fasta, "A")

#Load the mutational data for KRAS. Here the mutational data was obtained from the
#COSMIC database (version 58). 
data(KRAS.Mutations)

#Show the remapped order using the MDS remapper.
get.Remapped.Order(KRAS.Mutations, KRAS.Positions$Positions)

#Show the remapped order using the Linear remapper. Note that Amino Acid 61 is missing 
#if get.AlignedPositions was used. Thus Amino Acid 61 is missing in the returned result.
get.Remapped.Order(KRAS.Mutations, KRAS.Positions$Positions, method = "Linear")
												
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Order}
\keyword{Mutations}
\keyword{Amino Acids}

