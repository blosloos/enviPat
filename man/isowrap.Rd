\name{isowrap}
\alias{isowrap}

\title{
Combined (batch) calculation of isotope pattern, envelope and centroids/intensoids/valleys on
interpolated resolutions.
}
\description{
Wrapper combining the functions \code{\link[enviPat]{getR}}, \code{\link[enviPat]{isopattern}}, \code{\link[enviPat]{envelope}} and
\code{\link[enviPat]{vdetect}}.

Uses chemical formulas from \code{\link[enviPat]{check_chemform}} as argument.

}

\usage{
isowrap(isotopes, checked, resmass, resolution = FALSE, nknots = 6, 
spar = 0.2, threshold = 0.1, charge = 1, emass = 0.00054858, algo=2, 
ppm = FALSE, dmz = "get", frac = 1/4, env = "Gaussian", 
detect = "centroid", plotit = FALSE,verbose = TRUE )
}

\arguments{

  	\item{isotopes}{Dataframe listing all relevant isotopes, such as \code{\link[enviPat]{isotopes}}.}

  	\item{checked}{Output dataframe from \code{\link[enviPat]{check_chemform}} with correct chemical formulas.}

  	\item{resmass}{For resolution interpolation: dataframe with two columns, resolution and mass; see \code{\link[enviPat]{getR}}. 
			Otherwise, set to \code{FALSE} and use argument \code{resolution} to utilize a single resolution.} 

  	\item{resolution}{Single resolution value. Only used if argument \code{resmass} is set to \code{FALSE}.}

  	\item{nknots}{Number of knots, see \code{\link[enviPat]{getR}}. Ignored if argument \code{resmass} set to \code{FALSE}.}

  	\item{spar}{Smoothing parameter, see \code{\link[enviPat]{getR}}. Ignored if argument \code{resmass} set to \code{FALSE}.}

  	\item{threshold}{Abundance below which isotope peaks are omitted, see \code{\link[enviPat]{isopattern}}.}

  	\item{charge}{z in m/z, see \code{\link[enviPat]{isopattern}}.}

  	\item{emass}{Electrone mass. Only relevant if \code{charge} is not set to \code{FALSE}, see \code{\link[enviPat]{isopattern}}.}
	
  	\item{algo}{Which algorithm to use? Type \code{1} or \code{2}. See details section in \code{\link[enviPat]{isopattern}}.}

	\item{ppm}{Set stick discretization, see details section of \code{\link[enviPat]{envelope}}.}

  	\item{dmz}{Set stick discretization, see details section of \code{\link[enviPat]{envelope}}.}

  	\item{frac}{Set stick discretization, see details section of \code{\link[enviPat]{envelope}}.}

  	\item{env}{Peak shape function, see \code{\link[enviPat]{envelope}}.}

  	\item{detect}{Return either \code{"centroid"}, \code{"intensoid"} or \code{"valley"}. See \code{\link[enviPat]{vdetect}}.}

  	\item{plotit}{Should results be plotted, \code{TRUE/FALSE}?}
	
	\item{verbose}{Verbose, \code{TRUE/FALSE}?}
	
}



\value{
List with length equal to length of list \code{profiles}, with equal names of list entries.
Each entry in that list contains the centroids, intensoids or valley of the envelope in two columns:
\item{m/z}{m/z}
\item{abundance}{area(centroid) or abundance (intensoid, valley)}
}


\author{
Martin Loos, Christian Gerber
}


\seealso{
\code{\link[enviPat]{vdetect}}
}

\examples{

data(isotopes);
data(resolution_list);
data(chemforms);
chemforms<-chemforms[1:10];

checked<-check_chemform(
    isotopes,
    chemforms
  );

resmass<-resolution_list[[1]]

centro<-isowrap(
  isotopes,
  checked,
  resmass=resolution_list[[4]],
  resolution=FALSE,
  nknots=4,
  spar=0.2,
  threshold=0.1,
  charge=1,
  emass=0.00054858,
  algo=2,
  ppm=FALSE,
  dmz="get",   # retrieve dm from R=m/dm
  frac=1/4,
  env="Gaussian",
  detect="centroid",
  plotit=TRUE
)


}

