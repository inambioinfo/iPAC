\name{NEWS}
\title{NEWS file for the iPAC package}
\section{Release Notes for version 0.3.0 - 2012.06.26}{
  \itemize{
    \item First release of the iPAC package.
    \item ClusterFind method allows the use of MDS and linear mappers.
    \item Two beta methods available to reconcile data between the COSMIC
    and PDB databases.
  }
}

\section{Release Notes for version 1.1.2- 2012.10.14}{
  \itemize{
    \item Included get.Remapped.Order() function that displays the reshuffled amino acids after culling.
    \item Set the "AtomCount" column label to display "Can.Count" when the get.Positions function is called to be consistent with the get.AlignedPositions function.
	\item Included a "Plot.Protein.Linear" function that helps visualize the protein rearrangements.
	}
}

\section{Release Notes for version 1.1.2- 2012.10.17}{
  \itemize{
    \item You can now specify the title of your choice to the "Plot.Protein.Linear" function.
	}
}
\section{Release Notes for Version 1.1.3 - 2013.06.16}{
  \itemize{
    \item COSMIC changed it's display of FASTA sequences. Examples changed to work with the Uniprot Site.
  }
}

\section{Release Notes for Version 1.5.2 - 2013.06.26}{
  \itemize{
    \item Updated file to be consistent with Bioconductor version numbering. Also fixed a rare bug that gave an error if every residue aligned between mutational and positional data.
  }
}
