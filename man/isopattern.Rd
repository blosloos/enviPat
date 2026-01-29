\name{isopattern}

\alias{isopattern}

\title{Isotope pattern calculation}

\description{
The function calculates the isotopologues ("isotope fine structure") of a given chemical formula or a
set of chemical formulas (batch calculation) with fast and memory efficient transition tree algorithms,
which can handle relative pruning thresholds. 
Returns accurate masses, probabilities and isotopic compositions of individual isotopologues. 
The isotopes of elements can be defined by the user.
}

\usage{isopattern(isotopes, chemforms, threshold = 0.001, charge = FALSE, 
emass = 0.00054857990924, plotit = FALSE, algo=1, rel_to = 0, verbose = TRUE,
return_iso_calc_amount = FALSE)}

\arguments{
  
	\item{isotopes}{Dataframe listing all relevant isotopes, such as \code{\link[enviPat]{isotopes}}.}
  
	\item{chemforms}{Vector with character strings of chemical formulas, such as data set \code{\link[enviPat]{chemforms}} 
		or the second column in the value of \code{\link[enviPat]{check_chemform}}.}

  	\item{threshold}{Probability below which isotope peaks can be omitted, as specified by argument \code{rel_to}. 
		Set to \code{0} if all peaks shall be calculated.}

  	\item{charge}{z in m/z. Either a single integer or a vector of integers with length equal to that of argument \code{chemforms}. 
		Set to \code{FALSE} for omitting any charge calculations.}

  	\item{emass}{Electrone mass; only relevant if \code{charge} is not set to \code{FALSE}.}

  	\item{plotit}{Should results be plotted, \code{TRUE/FALSE}?}

	\item{algo}{Which algorithm to use? Type \code{1} or \code{2}. See details.}

	\item{rel_to}{Probability definition, numeric \code{0,1,2,3 or 4}? See details.}

	\item{verbose}{Verbose, \code{TRUE/FALSE}?}
	
	\item{return_iso_calc_amount}{Ignore; number of intermediate isotopologues.}
}



\details{

Isotope pattern calculation can be done by chosing one of two algorithms, set by argument \code{algo}. Both algorithms use
transition tree updates to derive the exact mass and probability of a new isotopologue from existing ones, by steps of single isotope replacements. 
These transition tree approaches are memory-efficient and fast for a wide range of molecular formulas and are able to reproduce the isotope fine structure
of molecules. The latter must often be pruned during calculation, c.p. argument \code{rel_to}. 

\code{algo=1} grows transition trees within element-wise sub-molecules, whereas \code{algo==2} grows them in larger sub-molecules of two elements, if available.
The latter approach can be slightly more efficient for very large or very complex molecules. The sub-isotopologues within sub-molecules are finally combined to
the isotopologuees of the full molecule. In contrast, intermediate counts of sub-isotopologues instead of fine structures are returned for \code{return_iso_calc_amount==TRUE}

\code{rel_to} offers 5 possibilities of how probabilities are defined and pruned, each affecting the \code{threshold} argument differently.
Default option \code{rel_to=0} prunes and returns probabilities relative to the most intense isotope peak; 
\code{threshold} states a percentage of the intensity of this latter peak.
Similarly, option \code{rel_to=1} normalizes relative to the peak consisting of the most abundant isotopes for each element, which
is often the monoisotopic one.
Option \code{rel_to=2} prunes and returns absolute probabilities ; \code{threshold} is not a percentage but an abolute cutoff.
Options \code{rel_to=3} and \code{rel_to=4} prune relative to the most intense and "monoisotopic" peak, respectively. 
Although \code{threshold} is a percentage, both options return absolut probabilities .

}


\value{

List with length equal to length of vector \code{chemforms}; names of entries in list = chemical formula in chemform.
Each entry in that list contains information on individual isotopologues (rows) with columns:

\item{m/z}{First column; m/z of an isotope peak.}
\item{abundance}{Second column; abundance of an isotope peak. Probabilities  are set relative to the most abundant peak of the isotope pattern.}
\item{12C, 13C, 1H, 2H, ...}{Third to all other columns; atom counts of individual isotopes for an isotope peak.}

}

\references{
Loos, M., Gerber, C., Corona, F., Hollender, J., Singer, H. (2015).  Accelerated isotope fine structure calculation using 
pruned transition trees, Analytical Chemistry 87(11), 5738-5744.

	\url{
		https://pubs.acs.org/doi/abs/10.1021/acs.analchem.5b00941
	}
	
	\url{
		https://www.envipat.eawag.ch/index.php
	}
}

\author{Martin Loos, Christian Gerber}

\note{
It is highly recommended to check argument \code{chemforms} with \code{\link[enviPat]{check_chemform}} prior to running 
\code{\link[enviPat]{isopattern}}; argument \code{chemforms} must conform to chemical formulas as defined in \code{\link[enviPat]{check_chemform}}.
Element names must be followed by numbers (atom counts of that element),  i.e. C1H4 is a valid argument whereas CH4 is not.
Otherwise, numbers may only be used in square brackets to denote individual isotopes defined in the element name column of iso_list, such as [14]C or [18]O.
For example, [13]C2C35H67N1O13 is the molecular formula of erythromycin labeled at two C-positions with [13]C;
C37H67N1O13 is the molecular formula of the unlabeled compound.

For correct adduct isotope pattern calculations, please check \code{\link[enviPat]{adducts}}.
}

\section{warning}{

Too low values for \code{threshold} may lead to unnecessary calculation of low probable isotope peaks - to the extent that not enough memory is available 
for either of the two algorithms.

}


\seealso{
	\code{\link[enviPat]{isopattern}}
	\code{\link[enviPat]{chemforms}}
	\code{\link[enviPat]{check_chemform}}
	\code{\link[enviPat]{getR}}
	\code{\link[enviPat]{envelope}}
	\code{\link[enviPat]{vdetect}}
	\code{\link[enviPat]{check_several}}
}

\examples{

############################
# batch of chemforms #######
data(isotopes)
data(chemforms)
pattern<-isopattern(
  isotopes,
  chemforms,
  threshold=0.1,
  plotit=TRUE,
  charge=FALSE,
  algo=1
)
############################
# Single chemical formula ##
data(isotopes) 
pattern<-isopattern(
  isotopes,
  "C100H200S2Cl5",
  threshold=0.1,
  plotit=TRUE,
  charge=FALSE,
  algo=1
)
############################


}



