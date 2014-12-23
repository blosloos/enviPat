\name{check_ded}
\alias{check_ded}

\title{
Check if a chemical formula is subset in another one
}
\description{
Check if a chemical formula is contained in another chemical formula
}
\usage{
check_ded(formulas, deduct)
}

\arguments{
  \item{formulas}{Vector with the containing chemical formula(s)}
  \item{deduct}{Chemical formula to be contained ("deduct")}
}

\value{
Returns a vector with length 0f input \code{formulas}, with \code{TRUE} if deduct is not contained and \code{FALSE} otherwise.
}

\author{Martin Loos}

\note{
Might be used used prior to subtracting a "deduct" chemical formula from that of a molecule when including adducts in the calculation of 
isotopic patterns. Chemical formulas must conform to what is described in \code{\link[enviPat]{check_chemform}}.
}


\seealso{\code{\link[enviPat]{adducts}}}

\examples{
formulas<-c("C8H4Cl2","C10H16O2","C3H10")
deduct<-c("C4H10")
check_ded(formulas, deduct)

}