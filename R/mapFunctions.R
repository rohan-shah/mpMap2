#' @title Map functions
#' @name mapFunctions
NULL
#' @describeIn mapFunctions Convert from Haldane distance to recombination fraction
#' @param x centiMorgan distance
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
#' @export
rfToHaldane <- function(r) 
{
	return(-50*log(1-2*r))
}
#' @describeIn mapFunctions Convert from recombination fraction to Kosambi distance
#' @param r recombination fraction
#' @export
rfToKosambi <- function(r) 
{
	return(50*atanh(2*r))
}
#' @describeIn mapFunctions Convert from Kosambi distance to recombination fraction
#' @param x centiMorgan distance
#' @export
kosambiToRf <- function(x)
{
	return(.5*tanh(2*x/100))
}
#' @describeIn mapFunctions Convert from recombination fraction to Kosambi distance
#' @export
kosambi <- kosambiToRf
