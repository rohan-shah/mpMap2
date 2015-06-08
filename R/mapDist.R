#' Conversion between recombination fractions (R) and map distance (X) 
#'
#' Haldane and Kosambi map functions and inverses for converting between 
#' recombination fractions and map distance (in cM)
#' @export
#' @rdname mapDist
#' @aliases haldaneX2R haldaneR2X kosambiR2X kosambiX2R
#' @param x Map distance (measured in centiMorgans)
#' @param r Recombination fraction
#' @references Jurg Ott, Analysis of Human Genetic Linkage, Johns Hopkins University Press, Baltimore 1999

haldaneR2X <-
function(r) return(-50*log(1-2*r))

#' @rdname mp-mapdist
haldaneX2R <-
function(x) return(.5*(1-exp(-2*x/100)))

#' @rdname mp-mapdist
kosambiR2X <-
function(r) return(50*atanh(2*r))

#' @rdname mp-mapdist
kosambiX2R <-
function(x) return(.5*tanh(2*x/100))