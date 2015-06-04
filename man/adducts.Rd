\name{adducts}
\alias{adducts}
\docType{data}
\title{
Adduct list
}
\description{
List of common adducts observed for ESI-MS measurements in soft positive and negative ionization modes.
}
\usage{data(adducts)}
\format{
  A data frame with 47 observations on the following 6 variables.
  \describe{
    	\item{\code{Name}}{Adduct name}
    	\item{\code{calc}}{Equation for calculating adduct m/z from uncharged non-adduct molecular mass M (m/z = M/z + X)}
    	\item{\code{Charge}}{z}
    	\item{\code{Mult}}{1/z}
    	\item{\code{Mass}}{X}
    	\item{\code{Ion_mode}}{Ionization mode (positive or negative)}
	\item{\code{Formula_add}}{Adduct chemical formula to be added}
	\item{\code{Formula_ded}}{Adduct chemical formula to be subtracted}
	\item{\code{Multi}}{Factor to multiply chemical formula with}
  }
}

\details{
The correct way to calculate the isotopic pattern of a specific adduct is the following.
First, multiply the chemical formula of the molecule by the times it appears in the final adduct; \code{\link[enviPat]{multiform}}.  
Second, add the chemical formula of any adduct to that of the molecule; \code{\link[enviPat]{mergeform}}.
Third, subtract the chemical formula of any deduct from that of the molecule; \code{\link[enviPat]{check_ded}} & \code{\link[enviPat]{subform}}.
Finally, calculate the isotopic fine structure using the correct \code{charge} argument in \code{\link[enviPat]{isopattern}}.   

}

\note{
Chemical formulas must conform to what is described in \code{\link[enviPat]{check_chemform}}.
}

\seealso{
	\code{\link[enviPat]{multiform}}
	\code{\link[enviPat]{mergeform}}
	\code{\link[enviPat]{check_ded}}
	\code{\link[enviPat]{subform}}
}

\source{
\url{http:/fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/}
}

\references{
Huang N., Siegel M.M., Kruppa G.H., Laukien F.H., J. Am. Soc. Mass. Spectrom. 1999, 10.
Automation of a Fourier transform ion cyclotron resonance mass spectrometer for acquisition, analysis, 
and e-mailing of high-resolution exact-mass electrospray ionization mass spectral data
}

\examples{
# example of M+H adduct batch calculation
data(adducts)
data(isotopes)
data(chemforms)
# (1) check formulas for consistency - recommended
checked_chemforms <-	check_chemform(isotopes, chemforms)			 
# (2) multiply, see column 4 of adducts
chemforms <-multiform(checked_chemforms[,2],1)			 
# (3) add adduct - see column 7 of adducts
chemforms<-mergeform(chemforms,"H1")
# (4) calculate fine structure
patterns <- isopattern(isotopes, chemforms)
}
\keyword{datasets}
