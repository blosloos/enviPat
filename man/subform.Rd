\name{subform}
\alias{subform}

\title{Subtract one chemical formula from another}

\description{Subtract one chemical formula from another}


\usage{subform(formula1,formula2)}

\arguments{

  \item{formula1}{Chemical formula to subtract from}

  \item{formula2}{Chemical formula to subtract}

}

\details{
Useful for adduct calculations, check \code{\link[enviPat]{adducts}}.
Chemical formulas must conform to what is described in \code{\link[enviPat]{check_chemform}}.
Prior check if \code{formula2} is contained in \code{formula2} at all? See \code{\link[enviPat]{check_ded}}.
}

\value{A unified and filtered peaklist}

\author{Martin Loos}

\seealso{
\code{\link[enviPat]{adducts}},\code{\link[enviPat]{check_ded}}
}


\examples{
formula1<-c("C10[13]C2H10Cl10")
formula2<-c("C2H5[13]C1")
subform(formula1,formula2)
}











