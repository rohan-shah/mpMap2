#' @export
orderCross <- function(mpcrossLG)
{
	if(!is(mpcrossLG, "mpcrossLG"))
	{
		stop("Input object must have linkage groups")
	}
	mpcrossLG <- as(mpcrossLG, "mpcrossLG")
	permutation <- .Call("order", mpcrossLG, mpcrossLG@lg@allGroups, PACKAGE="mpMap2")
	return(subset(mpcrossLG, markers = permutation))
}
