\name{nmc}
\alias{nmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
nmc
}
\description{
This runs the clustering algorithm developed by Ye. et. al on the original linear protein.
}
\usage{
nmc(x, alpha = 0.05, multtest = c("Bonferroni", "BH", "None"))
}
\arguments{
  \item{x}{
	x is a matrix of 0's and 1's. Each column represents a sample (cell-line, individual, etc.) while each column represents an amino acid in the protein. The i-th column represents the i-th amino acid.
}
  \item{alpha}{
	The significance level to be used. 
}
  \item{multtest}{
	The type of multiple comparisons adjustment to be used.
}
}
\value{
Returns a matrix of results which include the starting and ending positions of the clusters along with a p-value for each cluster.
}
\references{
Ye et. al., Statistical method on nonrandom clustering with application to somatic mutations in cancer. \emph{BMC Bioinformatics}. 2010. doi:10.1186/1471-2105-11-11.
}

\note{
Please see the original paper and methodology for a more detailed description of each of the parameters to this function.
}

\examples{
data(KRAS.Mutations)
nmc(KRAS.Mutations, alpha = 0.05, multtest = "Bonferroni")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{NMC}
