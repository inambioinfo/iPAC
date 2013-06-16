\name{get.Positions}
\alias{get.Positions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.Positions
}
\description{
This function reads a CIF file to extract the names and (x,y,z) coordinates of each residue. Then, the information within the CIF file itself is used to identify the canonical position for each amino acid.
}
\usage{
get.Positions(CIF.File.Location, Fasta.File.Location, chain.required = "A",
			  RequiredModelNum = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CIF.File.Location}{
		The location of the CIF file to be read.
}
  \item{Fasta.File.Location}{
	The location of the FASTA (or FASTA-like) file to be read.
}
  \item{chain.required}{
	The side chain in the protein from which to extract positions in the CIF file.
}
  \item{RequiredModelNum}{
	The required model num to extract positions from in the CIF file. If the RequiredModelNum == NULL, the method will use the first model number found in the file.
}
}
\details{
This algorithm attempts to recreate the canonical sequence by considering 3 seperate parts of a CIF file. The first part, the "_pdbx_poly_seq_scheme" section corresponds to the SEQRES entries in PDB files and lists all the amino acid residues in the protein that were studied including those that were not included in the 3D structure. The second part, the "struct_ref_seq" section, links the Uniprot (which provides the canonical sequence), PDB and internal CIF amino acid numbering schemes. Finally, the "_atom_site" section provides the (x,y,z) positional information. All three parts are tied together via a "seq_id" numerical identifier. It is thus possible to use the numbering scheme map provided in the "struct_ref_seq" section to calculate the canonical numbering for the residues listed in the "_pdbx_poly_seq_scheme" and then finally match them up to the residues listed in the "_atom_site" section. 
}
\value{
\item{Positions}{A dataframe that shows the amino acids extracted, their position in protein ordering, the side chain and their position in 3D space.}
\item{External.Mismatch}{This shows the mismatches that were found between the amino acids extracted from the file and the canonical sequence in the FASTA file.}
\item{PDB.Mismatch}{This shows the mismatches that are known to exist using only the information available in the CIF file.}
\item{Result}{The final status of extracting the amino acid sequence. It will either display "OK" or an error message.}
}
\note{
Due to the complexity of CIF file, this function is only able to account for a subset of the structures in the PDB database. If this function fails, consider using the the \emph{get.AlignedPositions} function. 

If you would like to contribute to this function, please contact the author.
}
\examples{
###################################################
#Observe that the final result from the code below is "OK". That is because the only
#mismatched residue at position 61, was documented in the CIF file as well.
#Thus it is considered a "reconciled" mismatch. It is up to the user to decide if 
#they want to include it in the position sequence or remove it. 

CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.extracted.positions<- get.Positions(CIF, Fasta, "A")
###################################################


###################################################
#Observe that the final result from the code below is "OK". Here we use a different file
#location for the canonical sequence -- the UNIPROT database. The canonical sequence is slightly 
#different and matches up exactly to the extracted positions. As there is only 1 isoform listed
#on UNIPROT for PIK3CA we suggest that users obtain the mutational data and the canonical 
#sequence information from the same source. For example, if your mutation data was obtained from 
#COSMIC, you should use COSMIC to get the canonical protein sequence.

CIF <- "http://www.pdb.org/pdb/files/2RD0.cif"
Fasta <- "http://www.uniprot.org/uniprot/P42336.fasta"
PIK3CA.extracted.positions<- get.Positions(CIF, Fasta, "A")
###################################################
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Amino Acids}
\keyword{Positions}
\keyword{CIF}
