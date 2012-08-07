\name{get.AASeq}
\alias{get.AASeq}
\title{
get.AASeq
}
\description{
Reads a FASTA (or similarly formatted) file and returns all the amino acids as a list.
}
\usage{
get.AASeq(filelocation)
}
\arguments{
  \item{filelocation}{
	The location of the file being read. 
}
}
\details{
The first line in the file being processed must begin with a ">". If the ending character in the file is a "*" (indicating a terminator -- as is common in the COSMIC database), this character is excluded from the sequence as well.
}
\value{
	An ordered list of amino acids in the protein. Each element of the list represents an amino acid.
}  
\examples{
#Reads the amino acid sequence for KRAS from COSMIC.
get.AASeq("http://www.sanger.ac.uk/perl/genetics/CGP/cosmic_sequence?ln=KRAS;format=aa")
}
\keyword{Amino Acids}
\keyword{FASTA}
