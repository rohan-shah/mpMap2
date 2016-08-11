#' @export
as.mpInterval <- function(object)
{
	if(!isS4(object) || !is(object, "mpcrossMapped"))
	{
		stop("Input object must be an S4 object of class mpcrossMapped")
	}
}
