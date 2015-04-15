\name{mergeform}
\alias{mergeform}

\title{Combine chemical formulas}

\description{Combine chemical formulas}


\usage{mergeform(formula1,formula2)}

\arguments{

  \item{formula1}{Vector of first chemical formula(s), character string(s)}

  \item{formula2}{Second chemical formula, single character string}

}

\details{
Useful for adduct calculations, check \code{\link[enviPat]{adducts}}.
Chemical formulas must conform to what is described in \code{\link[enviPat]{check_chemform}}.
}


\value{Merged chemical formula(s), character string}

\author{Martin Loos}

\seealso{
\code{\link[enviPat]{adducts}}
}


\examples{
formula1<-c("C10[13]C2H10Cl10")
formula2<-c("C2H5Na1")
mergeform(formula1,formula2)
}











