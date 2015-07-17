#' @title Map functions
#' @name mapFunctions
NULL
#' @describeIn mapFunctions Convert from Haldane distance to recombination fraction
#' @export
haldaneToRf <- function(x)
{
	return(.5*(1-exp(-2*x/100)))
}
#' @describeIn mapFunctions Convert from recombination fraction to Haldane distance
#' @export
rfToHaldane <- function(r) 
{
	return(-50*log(1-2*r))
}
#' @describeIn mapFunctions Convert from recombination fraction to Kosambi distance
#' @export
rfToKosambi <- function(r) 
{
	return(50*atanh(2*r))
}
#' @describeIn mapFunctions Convert from Kosambi distance to recombination fraction
#' @export
kosambiToRf <- function(x)
{
	return(.5*tanh(2*x/100))
}