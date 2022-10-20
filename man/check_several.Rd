\name{check_several}

\alias{check_several}

\title{Check for overlapping molecules.}

\description{
	Check for molecules overlapping in m/z, 
	based on isotope fine structures from \code{\link[enviPat]{isopattern}}
	or on centroids/intensoids from \code{\link[enviPat]{envelope}}.
}

\usage{check_several(pattern, dmz, ppm = TRUE)}

\arguments{

  \item{pattern}{Output from \code{\link[enviPat]{isopattern}} or from \code{\link[enviPat]{envelope}}.}

  \item{dmz}{m/z window. In combination with \code{ppm=TRUE} set as ppm or with \code{ppm=FALSE} set as absolute m/z.}

  \item{ppm}{Should m/z window be set in ppm (\code{TRUE}) or absolute m/z (\code{FALsE})?}

}

\details{
Overlaps in m/z among molecules are screened for within the m/z tolerance defined by the arguments \code{dmz} and \code{ppm}. 
}

\value{
Dataframe with 4 columns, with number of rows equal to the length of argument \code{pattern}

	\item{compound}{Chemical formula of the compound}
	\item{warning}{Overlap detected?}
	\item{to?}{If overlap: with wich other compound(s)? Refers to row number, recycled for peak_number.}
	\item{peak_number}{If overlap: with which peak(s) of the other compound(s)? Refers to peak number.}

}

\author{Martin Loos, Christian Gerber}


\seealso{
	\code{\link[enviPat]{isopattern}}
	\code{\link[enviPat]{envelope}}
}

\examples{
data(isotopes)
data(chemforms)
pattern<-isopattern(
  isotopes,
  chemforms,
  threshold=0.1,
  plotit=TRUE,
  charge=FALSE,
  emass=0.00054858,
  algo=1
)
check_several(pattern,dmz=0.001,ppm=FALSE)
}

