#' @title Map functions
#' @name mapFunctions
NULL
#' @describeIn mapFunctions Convert from Haldane distance to recombination fraction
#' @description Functions used to convert between recombination fractions and centiMorgan distances. 
#' @param x centiMorgan distance
#' @return Recombination fraction.
#' @export
haldaneToRf <- function(x)
{
	return(.5*(1-exp(-2*x/100)))
}
#' @describeIn mapFunctions Convert from Haldane distance to recombination fraction
#' @export
haldane <- haldaneToRf
#' @describeIn mapFunctions Convert from recombination fraction to Haldane distance
#' @param r recombination fraction
#' @return Genetic distance in cM. 
#' @export
rfToHaldane <- function(r) 
{
	return(-50*log(1-2*r))
}
#' @describeIn mapFunctions Convert from recombination fraction to Kosambi distance
#' @return Genetic distance in cM. 
#' @export
rfToKosambi <- function(r) 
{
	return(50*atanh(2*r))
}
#' @describeIn mapFunctions Convert from Kosambi distance to recombination fraction
#' @return Recombination fraction.
#' @export
kosambiToRf <- function(x)
{
	return(.5*tanh(2*x/100))
}
#' @describeIn mapFunctions Convert from recombination fraction to Kosambi distance
#' @return Recombination fraction.
#' @export
kosambi <- kosambiToRf
