\name{check_chemform}

\alias{check_chemform}

\title{Check chemical formulas}

\description{Checks chemical formulas (=a vector of character strings) for consistency with usage in \code{\link[enviPat]{isopattern}}; 
 calculation of the molecular mass.}

\usage{check_chemform(isotopes,chemforms,get_sorted=FALSE)}

\arguments{
  \item{isotopes}{\code{\link[enviPat]{isotopes}}}
  \item{chemforms}{Vector of character strings with chemical formulas}
  \item{get_sorted}{Should elements in each formula be sorted according to their order in \code{isotopes}?}
}


\details{

	Checks if (1) a chemical formula contains only letters, numbers and square brackets, (2) elements can be found 
	in \code{\link[enviPat]{isotopes}} and (3) letters are all followed by a number of atom counts.
	If (3) is missing, it is automatically set to 1.
	

	(2) must consist of an upper case letter, possibly followed by lower case letters; to refer to individual isotopes
	(e.g., from isotope labelling of a molecule, e.g., N5 vs. [15]N2N3), square brackets may precede the capital letter.
	Anything else that may usually be part of a standard chemical formula (e.g., (+), (C4H3)2, dashes,...) is not permissible.


	The molecular mass will be calculated from isotope masses and abundances listed in \code{\link[enviPat]{isotopes}}.

}


\value{
Dataframe with 3 columns:
\item{warning}{Correct chemical formula, \code{FALSE}/\code{TRUE}?}
\item{new_formula}{Chemical formula}
\item{monoisotopic_mass}{Monoisotopic mass}
}


\author{
Martin Loos, Christian Gerber
}


\note{
Highly recommended for usage with \code{\link[enviPat]{isopattern}}
}


\seealso{
\code{\link[enviPat]{isopattern}}
\code{\link[enviPat]{isotopes}}
}


\examples{

# Check package data set of chemical formulas #############
data(chemforms);
data(isotopes);
checked<-check_chemform(isotopes,chemforms);
checked;

# Check for some senseless molecular formulas #############
chemforms<-c("C900Cl4H49","O82394","C8G500Zn9","Br1","6DBr9889");
data(isotopes);
checked<-check_chemform(isotopes,chemforms);
checked;

# Molecular mass with and without isotope labelling #######
chemforms<-c("C10H5N4O5","[13]C2C8D2H3[15]N2N2[18]O2O3");
data(isotopes);
checked<-check_chemform(isotopes,chemforms);
checked;

}




















