\name{enviPat-package}
\alias{enviPat-package}
\alias{enviPat}
\docType{package}
\title{ Calculation of isotope patterns, stick profiles (envelopes) and centroids/intensoids for mass spectrometry.}

\description{

Fast and memory-efficient calculation of isotope patterns (fine structures) for up to very large molecules, based on three different algorithms.
Subsequent convolution of isotope patterns with a peak shape function to theoretical envelopes (profiles).
Based on envelopes, valley detection and centroidization/intensoid calculation.
Allows for batch processing of chemical formulas and interpolation of measurement resolutions.
Includes a wrapper combining all of the above functionalities.


Furthermore, includes 
(1) a check for consistency of chemical formulas,
(2) a check for molecules with overlapping isotope patterns,
(3) a list of all stable isotopes,
(4) a list of different resolution data sets for Thermo Orbitrap and QExactive high-resolution mass spectrometers and
(5) a list of adducts formed during electorspray ionization (ESI).


A web-based GUI for enviPat is freely available under \url{www.envipat.eawag.ch/}.
}

\details{
\tabular{ll}{
Package: \tab enviPat\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-05\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

}


\author{
	Martin Loos, Christian Gerber

	Maintainer: Martin Loos <loosmart@eawag.ch>
}

\keyword{ package }
\seealso{

	\code{\link[enviPat]{check_chemform}}
  	\code{\link[enviPat]{getR}}
  	\code{\link[enviPat]{isopattern}}
  	\code{\link[enviPat]{envelope}}
  	\code{\link[enviPat]{vdetect}}
  	\code{\link[enviPat]{isowrap}}
  	\code{\link[enviPat]{check_several}}

  	\code{\link[enviPat]{isotopes}}
  	\code{\link[enviPat]{resolution_list}}
  	\code{\link[enviPat]{chemforms}}

  	\code{\link[enviPat]{adducts}}
	\code{\link[enviPat]{check_ded}}
	\code{\link[enviPat]{mergeform}}
	\code{\link[enviPat]{subform}}
	\code{\link[enviPat]{multiform}}
}

