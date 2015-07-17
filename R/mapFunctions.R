#' @export
haldaneToRf <- function(x)
{
	return(.5*(1-exp(-2*x/100)))
}
#' @export
rfToHaldane <- function(r) 
{
	return(-50*log(1-2*r))
}
#' @export
rfToKosambi <- function(r) 
{
	return(50*atanh(2*r))
}
#' @export
kosambiToRf <- function(x)
{
	return(.5*tanh(2*x/100))
}