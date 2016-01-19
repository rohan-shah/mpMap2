#' @export
orderCross <- function(mpcrossLG, cool = 0.5, tmin = 0.1, nReps = 1)
{
	if(!is(mpcrossLG, "mpcrossLG"))
	{
		stop("Input object must have linkage groups")
	}
	mpcrossLG <- as(mpcrossLG, "mpcrossLG")
	permutation <- .Call("order", mpcrossLG, mpcrossLG@lg@allGroups, cool, tmin, nReps, PACKAGE="mpMap2")
	return(subset(mpcrossLG, markers = permutation))
}
