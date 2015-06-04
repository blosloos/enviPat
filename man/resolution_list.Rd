\name{resolution_list}
\alias{resolution_list}
\docType{data}
\title{Resolutions (R) list for Themo Orbitrap and QExactive mass spectrometers}
\description{

	List of different resolutions R=f(m/z) for various high-resolution mass spectrometers.
	For each of the instruments, different resolution settings are available. 
	Here, R is defined as R=(m/z)/(dm/z), with dm/z = peak width at half maximum.
	Serves as input to \link[enviPat]{getR} to interpolate R from given molecular masses.

}
\usage{data(resolution_list)}
\format{

The format is:
List with 19 data sets: Instrument/Resolution@m/z

Elite/R240000@400

Elite/R120000@400

Elite/R60000@400

Elite/R30000@400

OrbitrapXL,Velos,VelosPro/R120000@400

OrbitrapXL,Velos,VelosPro/R60000@400

OrbitrapXL,Velos,VelosPro/R30000@400

OrbitrapXL,Velos,VelosPro/R15000@400

OrbitrapXL,Velos,VelosPro/R7500@400

Q-Exactive,ExactivePlus/280K@200

Q-Exactive,ExactivePlus/R140000@200

Q-Exactive,ExactivePlus/R70000@200

Q-Exactive,ExactivePlus/R35000@200

Q-Exactive,ExactivePlus/R17500@200

Exactive/R100000@200

Exactive/R50000@200

Exactive/R25000@200

Exactive/R12500@200

OTFusion,QExactiveHF/480000@200

OTFusion,QExactiveHF/240000@200

OTFusion,QExactiveHF/120000@200

OTFusion,QExactiveHF/60000@200

OTFusion,QExactiveHF/30000@200

OTFusion,QExactiveHF/15000@200

TripleTOF5600/R28000@200 

QTOF_XevoG2-S/R25000@200
  
}
\source{
Data assembled from individual measurements.
}
\examples{
data(resolution_list)
resolution_list
}
\keyword{datasets}
