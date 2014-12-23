\name{getR}

\alias{getR}

\title{Interpolation of MS measurement resolution}

\description{
	Given a set of MS measurement resolutions (R) as a function of measurement mass (m/z), \link[enviPat]{getR} interpolates R for any given molecular mass(es)
	calculated by \code{\link[enviPat]{check_chemform}} using \code{\link[stats]{smooth.spline}}.
}

\usage{getR(checked, resmass, nknots = 13, spar = 0.1, plotit = TRUE)}

\arguments{

  \item{checked}{Dataframe produced by \code{\link[enviPat]{check_chemform}}.}

  \item{resmass}{Dataframe with two columns, resolution and mass; such as the list entries in \code{\link{resolution_list}}.}

  \item{nknots}{Integer number of knots to use for the smoothing spline. Default = 6. See also \code{\link[stats]{smooth.spline}}.}

  \item{spar}{Smoothing parameter, (0,1]. See also \code{\link[stats]{smooth.spline}}.}

  \item{plotit}{Plot results, \code{TRUE/FALSE} ?}

}

\value{Vector with resolutions.}

\note{
\code{\link[enviPat]{check_chemform}} gives molecular masses (m/z) for z=+/-1 only. 
If z>1 or z<-1 is required, molecular mass entries in argument checked have to be divided accordingly to be consistent.
}

\author{
Martin Loos, Christian Gerber
}


\seealso{
\code{\link[stats]{smooth.spline}}
\code{\link[enviPat]{check_chemform}}
\code{\link[enviPat]{resolution_list}}
}

\examples{

data(resolution_list)
resmass<-resolution_list[[4]]
data(isotopes)
data(chemforms)
checked<-check_chemform(isotopes,chemforms)
resolution<-getR(checked,resmass,nknots=13,spar=0.1,plotit=TRUE)

# same for z=-2:
checked<-check_chemform(isotopes,chemforms)
checked[,3]<-(checked[,3]/abs(-2))
resolution<-getR(checked,resmass,nknots=13,spar=0.1,plotit=TRUE)
}

