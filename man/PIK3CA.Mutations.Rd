\name{PIK3CA.Mutations}
\alias{PIK3CA.Mutations}
\docType{data}
\title{
PIK3CA.Mutations
}
\description{
 A matrix of the verified somatic mutations for PIK3CA from COSMIC.
}
\usage{data(PIK3CA.Mutations)}
\format{
  There are 96 rows, each representing a sample mutation and 1068 colums -- one for each amino acid in the protein.
}

\source{
Version 58 of the cosmic database was used.
}
\references{
S A Forbes, G Bhamra, S Bamford, E Dawson, C Kok, J Clements, A Menzies, J W Teague,
P A Futreal, and M R Stratton. The catalogue of somatic mutations in cancer (COSMIC).
Current Protocols in Human Genetics / Editorial Board, Jonathan L. Haines ... [et Al.],
Chapter 10:Unit 10.11, April 2008. ISSN 1934-8258. doi: 10.1002/0471142905.hg1011s57.

URL http://www.ncbi.nlm.nih.gov/pubmed/18428421. PMID: 18428421.
}
\examples{
data(PIK3CA.Mutations)

#Show the first 10 samples
PIK3CA.Mutations[1:10,]

#Show which samples had position 1047 mutated
PIK3CA.Mutations[,1047]
}
\keyword{datasets}
