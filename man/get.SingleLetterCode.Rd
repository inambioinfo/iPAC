\name{get.SingleLetterCode}
\alias{get.SingleLetterCode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.SingleLetterCode
}
\description{
This function takes in the 3 letter amino acid representation and returns the 1 letter representation. 
}
\usage{
get.SingleLetterCode(code)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{code}{
The 3 letter code in all capital letters.
}
}
\details{
If the 3 letter code is not found (either because it does not exist or because it is not in all capital letters), an "X" is returned.

This function is used as a helper function of \emph{get.AlignedPositions} and \emph{get.Positions}. It is not necessary to invoke it independently and we include it here solely as a matter of convenience for the user.
}
\value{
An amino acid single letter code returned as a string. The string will be an element in:

\{"P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T","X"\}.
}
\examples{
get.SingleLetterCode("MET")

get.SingleLetterCode("GLN")
 
get.SingleLetterCode("gln")

get.SingleLetterCode("DMY")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Amino Acids}