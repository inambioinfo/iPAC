\name{ClusterFind}
\alias{ClusterFind}
\title{
ClusterFind
}
\description{
ClusterFind is the main method of the \emph{iPAC} package. It identifies clusters of mutated amino acids while taking into account the protein structure.
}
\usage{
ClusterFind(mutation.data, position.data, method = "MDS", alpha = 0.05, 
			MultComp = "Bonferroni", Include.Culled = "Y", Include.Full = "Y", 
			create.map = "Y", Show.Graph = "Y", Graph.Output.Path = NULL,
			Graph.File.Name = "Map.pdf", Graph.Title = "Mapping", 
			OriginX = min(position.data[, 4]), OriginY = min(position.data[, 5]),
			OriginZ = min(position.data[, 6]))
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
  \item{alpha}{
	The significance level used in the NMC calculation. Please see \emph{Ye. et. al.} for more information.
}
  \item{MultComp}{
	The multiple comparisons adjustment used in the NMC calculation. Possible options are "None", "Bonferroni" and "BH". Please see \emph{Ye. et. al.} for more information.
}
  \item{Include.Culled}{
	If "Y", the standard NMC algorithm will be run on the protein after removing the amino acids for which there is no positional data.
}
  \item{Include.Full}{
	If "Y", the standard NMC algorithm will be run on the full protein sequence.
}
  \item{create.map}{
	If "Y", a graphical representation of the the dimension reduction from 3D to 1D space will be created (though not necessarily displayed).
}
  \item{Show.Graph}{
	If "Y", the graph representation will be displayed. Warning: You must be running R in a GUI environment, otherwise, an error will occur. 
}
  \item{Graph.Output.Path}{
	If you would like the picture saved atomatically to the disk, specify the output directory here. The \emph{Graph.File.Name} variable must be set as well.
}
  \item{Graph.File.Name}{
	If you would like the picture saved automatically to the disk, specify the output file name. The \emph{Graph.Output.Path} variable must be set as well.
}
  \item{Graph.Title}{
	The title of the graph to be created.
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
 The linear method fixes a point, defined by the parameters \emph{OriginX}, \emph{OriginY}, \emph{OriginZ}, and then calculates the distance from each amino acid to that point. The graph produced by \emph{ClusterFind} (if requested), shows these distances as dotted green lines. The length of the green line is used to reorder the protein, with the amino acid that corresponds to the shortest green line being ordered first and the amino acid corresponding to the longest green line being ordered last. 
 
 Additional methods will be available in future versions of this package. 
}
\value{
\item{Remapped}{This shows the clusters found while taking the 3D structure into account.}
\item{OriginalCulled}{This shows the clusters found if you run the NMC algorithm on the canonical linear protein, but with the amino acids for which we don't have 3D positional data removed.}
\item{Original}{This shows the clusters found if you run the NMC algorithn on the canonical linear protein with all the amino acids.}
\item{MissingPositions}{This shows which amino acids are present in the mutation matrix but for which we do not have positions. These amino acids are cut from the protein when calculating the \emph{Remapped} and \emph{OriginalCulled} results.}
}
\references{
Ye et. al., Statistical method on nonrandom clustering with application to somatic mutations in cancer. \emph{BMC Bioinformatics}. 2010. doi:10.1186/1471-2105-11-11.
}
\note{
If no significant clusters are found, a "NULL" will be returned for the appropriate section (Remapped, OriginalCulled, or Original).

If you want the graph to just display on a new R graphics devices without saving it to the disk, simply set the \emph{Graph.Output.Path} or the \emph{Graph.File.Name} parameters to be NULL while leaving both the \emph{create.map} and the \emph{Show.Graph} parameters to be "Y". 

If you are running this algorithm on a terminal with no GUI (such as a computational cluster), set \emph{Show.Graph} to "N" as R will not be able to open a new graphics device. However, you can still have the graphics saved for later viewing by setting the \emph{Graph.Output.Path} and \emph{Graph.File.Name} variables. 

If \emph{ClusterFind} displays the message "Error in 0:(n - j) : NA/NaN argument", this most likely signifies that after removing the amino acids for which there is no positional data, the mutation data matrix is all 0's. For instance, if all the mutations occured on the 5th amino acid in the protein, and we did not have 3D positional information for that amino acid, the mutation data for the remaining positions would be all 0's. An error message in this situation will be displayed in the results section as well. In such cases, the user should run the original NMC algorithm only (see \emph{nmc}) or select an alternative protein structure.

When unmapping back to the original space, the end points of the cluster in the mapped space are used as the endpoints of the cluster in the unmapped space.
}
\examples{
#Extract the data from a CIF file and match it up with the canonical protein sequence.
#Here we use the 3GFT structure from the PDB, which corresponds to the KRAS protein.
CIF<-"http://www.pdb.org/pdb/files/3GFT.cif"
Fasta<-"http://www.uniprot.org/uniprot/P01116-2.fasta"
KRAS.Positions<-get.Positions(CIF,Fasta, "A")

#Load the mutational data for KRAS. Here the mutational data was obtained from the
#COSMIC database (version 58). 
data(KRAS.Mutations)

#Identify and report the clusters using the default MDS method.
ClusterFind(mutation.data=KRAS.Mutations, 
							position.data=KRAS.Positions$Positions,
							create.map = "Y",Show.Graph = "Y")
							
#Identify and report the clusters using the linear method.	
ClusterFind(mutation.data=KRAS.Mutations, 
							position.data=KRAS.Positions$Positions,
							create.map = "Y",Show.Graph = "Y", method = "Linear")							
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Clusters}
\keyword{Mutations}
\keyword{Amino Acids}

