\name{Plot.Protein.Linear}
\alias{Plot.Protein.Linear}
\title{
Plot.Protein.Linear
}
\description{
Plot.Protein.Linear creates a visual representation of the protein reordering by showing the protein in linear form and setting the color based upon the remapped amino acid position.
}
\usage{
Plot.Protein.Linear(path, colCount, cex=0.5, height = 1,width = 1, color.palette = "heat")
}
\arguments{
  \item{path}{
	The remapped protein using only those amino acids for which mutational data is available.
}
  \item{colCount}{
	How many columns you want in the plot
}
  \item{cex}{
	How large you want the text to be. See R graphical parameters for more information.
}
  \item{height}{
	The height of the resulting box that outlines the amino acid number.
}
  \item{width}{
	The width of the resulting box that outlines the amino acid number.
}
  \item{color.palette}{
	One of of the standard color palettes available in R. The options are  "heat", "gray", "topo", and "cm".
}
}
\details{
  The terminal you are running this on must be able to create an R graphical device. 
}

\note{
This method is still in beta. If any bugs found, please email the package authors.
}

\examples{
#Extract the data from a CIF file and match it up with the canonical protein sequence.
#Here we use the 3GFT structure from the PDB, which corresponds to the KRAS protein.
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.sanger.ac.uk/perl/genetics/CGP/cosmic_sequence?ln=KRAS;format=aa"
KRAS.Positions<-get.Positions(CIF,Fasta, "A")

#Load the mutational data for KRAS. Here the mutational data was obtained from the
#COSMIC database (version 58). 
data(KRAS.Mutations)

#Show the remapped order using the MDS remapper.
new.ordering<-get.Remapped.Order(KRAS.Mutations, KRAS.Positions$Positions)

#Create the Plots. Note that Amino Acid 61 is missing if get.AlignedPositions
#was used. Thus Amino Acid 61 is missing in the plots below.
Plot.Protein.Linear(new.ordering, 25, color.palette = "heat")
Plot.Protein.Linear(new.ordering, 25, color.palette = "gray")
												
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Order}
\keyword{Mutations}
\keyword{Amino Acids}

