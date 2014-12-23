\name{multiform}
\alias{multiform}

\title{Multiply a chemical formula}

\description{Multiply all atom numbers in a chemical formula by a factor}


\usage{multiform(formula1,fact)}

\arguments{

  \item{formula1}{Chemical formula to be multiplied, vector of character strings}

  \item{fact}{Factor to multiply with}

}

\details{
Useful for adduct calculations, check \code{\link[enviPat]{adducts}}.
Chemical formulas must conform to what is described in \code{\link[enviPat]{check_chemform}}.
}

\value{Multiplied chemical formula, character string}

\author{Martin Loos}

\seealso{
\code{\link[enviPat]{adducts}}
}


\examples{
formula1<-c("C10[13]C2H10Cl10")
multiform(formula1,3)
}











