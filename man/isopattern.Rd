\name{isopattern}

\alias{isopattern}

\title{Isotope pattern calculation}

\description{
The function calculates the isotopic pattern (fine structure) of a given chemical formula or a
set of chemical formulas (batch calculation), based on three fast and memory efficient algorithms. 
The first algorithm can handle very large molecules and combinations of elements having many isotopes. 
Returns accurate masses, abundances and isotopic compositions of the individual isotopologues. 
The isotopes of elements can be defined by the user.
}

\usage{isopattern(isotopes, chemforms, threshold = 0.001, charge = FALSE, 
emass = 0.00054858, plotit = FALSE, algo=2, rel_to_mono = FALSE, verbose = TRUE)}

\arguments{
  
	\item{isotopes}{Dataframe listing all relevant isotopes, such as \code{\link[enviPat]{isotopes}}.}
  
	\item{chemforms}{Vector with character strings of chemical formulas, such as data set \code{\link[enviPat]{chemforms}} 
		or the second column in the value of \code{\link[enviPat]{check_chemform}}.}

  	\item{threshold}{Abundance below which isotope peaks can be omitted, given as percentage of the most abundant isotope peak of the molecule. 
		Set to \code{0} if all peaks shall be calculated.}

  	\item{charge}{z in m/z. Either a single integer or a vector of integers with length equal to that of argument \code{chemforms}. 
		Set to \code{FALSE} for omitting any charge calculations.}

  	\item{emass}{Electrone mass; only relevant if \code{charge} is not set to \code{FALSE}.}

  	\item{plotit}{Should results be plotted, \code{TRUE/FALSE}?}

	\item{algo}{Which algorithm to use? Type \code{1} or \code{2}. See details.}

	\item{rel_to_mono}{Should abundances be normalized relative to the monoisotopic instead of the most abundant peak, \code{TRUE/FALSE}?}

	\item{verbose}{Verbose, \code{TRUE/FALSE}?}
}



\details{

Isotope pattern calculation can be done by chosing one of three algorithms, set by argument \code{algo}. All algorithms use
hierarchical updates to derive the mass and abundance of a new isotopologue from an existing one, by steps of single isotope replacements.
Memory usage is lower and in most cases faster for the first two algorithms as compared to the third, 
allowing for calculation of very large molecules or inclusion of elements with many isotopes. Comparable in memory allocation, the second algorithm is 
faster for very small molecules than the first - but much slower for larger ones. 

The first algorithm \code{algo=1} uses tree-like combinatorial transitions to calculate daughter isotopologues from their parent node isotopologues, with
the monoisotopic composition as root node.
This approach first searches for branches of increasing abundance to find the isotopologue of maximum abundance, with transitions ordered as to minimize 
the occurrence of decreasing branches. The remaining branches are subsequently omitted if they (a) fall below a threshold relative to this most abundant 
isotopologue and (b) only contain branches of decreasing abundance. Furthermore, to avoid redundant calculations for transitions of the same isotope 
(but not the same isotopologues!), this global search is conducted in elementwise subtrees that are then combined.

The second algorithm \code{algo=2} does not use elementwise subtree maximum abundance searches with no differences to the first algorithm otherwise.

The third algorithm \code{algo=3} is similar to the one poposed by Li et al. (2010). Herein, mass states and abundances are calculated 
individually within seperate blocks for each of the elements present in a molecule without (!) abundance thresholds. 
These building-blocks are then combined to individual isotopologues, with peaks below the threshold abundance evenutally omitted. 
In the presented version, a fast calculation of elementwise building-blocks and their combination to isotopologues is implemented so as to avoid redundant 
calculations from both different updates or different combinations leading to the same isotopologue.

Note that when \code{rel_to_mono} is set to \code{TRUE}, the abundance \code{threshold} is specified relative to the monoisotopic instead of the most abundant peak.
}


\value{

List with length equal to length of vector \code{chemforms}; names of entries in list = chemical formula in chemform.
Each entry in that list contains information on individual isotopologues (rows) with columns:

\item{m/z}{First column; m/z of an isotope peak.}
\item{abundance}{Second column; abundance of an isotope peak. Abundances are set relative to the most abundant peak of the isotope pattern.}
\item{12C, 13C, 1H, 2H, ...}{Third to all other columns; atom counts of individual isotopes for an isotope peak.}

}

\references{

Loos, M. & Gerber, C., 201X. Tree-like hierarchy for the calculation of very large isotope patterns. To be submitted.

Li, L., Karabacak, N., Cobb, J., Wang, Q., Hong, P., agar, J., 2010. Memory-efficient calculation of the isotopic mass states of a molecule.
Rapid Communications in Mass Spectrometry, 24: 2689-2696.

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

Too low values for \code{threshold} may lead to unnecessary calculation of low abundance peaks - to the extent that not enough memory is available 
for either of the two algorithms. This is especially critial if \code{rel_to_mono} is set to \code{TRUE}.

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
  emass=0.00054858,
  algo=2
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
  emass=0.00054858,
  algo=2
)
############################


}



