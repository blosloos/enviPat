\name{isotopes}
\alias{isotopes}
\docType{data}
\title{
Stable isotopes
}
\description{
Dataframe with stable isotopes.
}
\usage{data(isotopes)}
\format{
  A data frame with 302 observations on the following 4 variables.
  \describe{
    \item{\code{element}}{Chemical element} 
    \item{\code{isotope}}{Stable isotopes of an element}
    \item{\code{mass}}{Relative atomic mass}
    \item{\code{abundance}}{Isotopic composition of an element}
	\item{\code{ratioC}}{Maximum number of atoms of an element for one C-atom in a molecule, based on 99.99 \% of case molecules.}
  }
}
\details{
The \code{ratioC}-value stems from a database survey conducted by Kind&Fiehn (2007); to disable, set value to 0.
The list serves as input into several package nontarget-functions. The first column of the data frame also 
contains names of specific isotopes used for labeled compounds. 
}
\source{
\url{http:/physics.nist.gov/cgi-bin/Compositions/stand_alone.pl}
}
\references{
Kind, T. and Fiehn, O., 2007. Seven golden rules for heuristic filtering  of molecular formulas obtained by 
accurate mass spectrometry. BMC Bioinformatics, 8:105.
}
\examples{
data(isotopes)
}
\keyword{datasets}
