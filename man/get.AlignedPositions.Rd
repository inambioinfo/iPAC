\name{get.AlignedPositions}
\alias{get.AlignedPositions}
\title{
get.AlignedPositions
}
\description{
This function reads a CIF file to extract the names and (x,y,z) coordinates of each residue. It then performs a pairwise alignment to convert the amino acid ordering in the CIF file to the canonical ordering specified by the FASTA file. The first element in the returned list, $Positions, is the positions matrix required by the ClusterFind method.
}
\usage{
get.AlignedPositions(CIF.File.Location, Fasta.File.Location, chain.required = "A",
		     RequiredModelNum = NULL, patternQuality = PhredQuality(22L),
		     subjectQuality = PhredQuality(22L), type = "global-local",
		     substitutionMatrix = NULL, fuzzyMatrix = NULL, gapOpening = -10,
		     gapExtension = -4, scoreOnly = FALSE)
}
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
  \item{patternQuality}{
    The patternQuality parameter in the pairwiseAlignment function in the Biostrings package.
}
  \item{subjectQuality}{
The subjectQuality parameter in the pairwiseAlignment function in the Biostrings package.
}
  \item{type}{
The type parameter in the pairwiseAlignment function in the Biostrings package. This should NOT be changed from "global-local" as we use the canonical protein from the FASTA file as the global pattern and the extracted positions from the CIF as the subject pattern. We then attempt to align parts of the subject pattern to the entire global pattern.
}
  \item{substitutionMatrix}{
The substitutionMatrix parameter in the pairwiseAlignment function in the Biostrings package.
}
  \item{fuzzyMatrix}{
The fuzzyMatrix parameter in the pairwiseAlignment function in the Biostrings package.
}
  \item{gapOpening}{
The gapOpening parameter in the pairwiseAlignment function in the Biostrings package.
}
  \item{gapExtension}{
The gapExtension parameter in the pairwiseAlignment function in the Biostrings package.
}
  \item{scoreOnly}{
The scoreOnly parameter in the pairwiseAlignment function in the Biostrings package.
}
}
\details{
This method is currently in BETA and is provided only as a convenient way to extract the required 3D positional information from a CIF file. Currently, CIF (and PDB) files can have a number of deviations from the canonical protein sequence including additional, missing and mismatched amino acids. The amino acid numbering in CIF files can also be different from the canonical protein. This makes it difficult to match up mutational data (from sources such as COSMIC) to 3D positional data (from sources such as the PDB) and necessitates the use of this function.

This method extracts the canonical amino acid sequence from the file at \emph{Fasta.File.Location}. It then attempts to align the amino acids extracted from the CIF file to the canonical sequence using the pairwiseAlignment function in the package \emph{Biostrings} that is available on \emph{Bioconductor}. After alignment, any amino acids that are mismatched between the canonical sequence and the extracted sequence are automatically removed so that the \emph{ClusterFind} method, which requires positional data as input, is only run on those amino acids which are correctly matched.
}
\value{
\item{Positions}{A dataframe that shows the extracted amino acids, their numerical position in the protein order, the protein side chain being used and the amino acid positions in 3D space.}
\item{Diff.Count}{A check that the amino acids remaining are in fact matched to the canonical protein. This returns a count of the number of amino acids remaining that do not match the canonical sequence and should be 0 if a successful alignment occurred.}
\item{Diff.Positions}{A description of the mismatched amino acids if any were found. If a succesful alignment occurred, this will be NULL.}
\item{Alignment.Result}{The raw alignment result returned by the pairwiseAlignment method in the \emph{Biostrings} package.}
\item{Result}{The final status of the alignment. If it is "OK", that means that this alignment appears to be ok. If the alignment failed, this item will contain an error message.}
}
\references{
\emph{Biostrings} Package. Bioconductor: Open software development for computational biology and bioinformatics R. Gentleman, V. J. Carey, D. M. Bates, B.Bolstad, M. Dettling, S.
  Dudoit, B. Ellis, L. Gautier, Y. Ge, and others 2004, Genome Biology, Vol. 5, R80
}
\note{
This function is provided as a convenience to the user so that they may easily obtain the required input for the \emph{ClusterFind} function. The user should validate the results.

The "Diff.Count" and "Diff.Positions" parts of the returned value should always be "0" and "NULL" if the positions were successfully extracted as mismatched positions are dropped automatically.  This behavior is different from the \emph{get.Positions} function where a reconciliation of differences is attempted using the information within the CIF file.

If you would like to contribute to this function, please contact the author.
}
\examples{
#Observe that position 61 is missing. It is atuomatically dropped as the pdb data
#specifies it as a "H" while the FASTA sequence specifies it as "Q".
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.sanger.ac.uk/perl/genetics/CGP/cosmic_sequence?ln=KRAS;format=aa"
get.AlignedPositions(CIF,Fasta, "A")
}
\keyword{Amino Acids}
\keyword{Positions}
\keyword{CIF}
