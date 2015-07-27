checkHets <- function(object)
{
	if(is.null(names(object)) || any(names(object) == ""))
	{
		return("Every entry must have a valid name")
	}
	if(length(unique(names(object))) != length(object))
	{
		return("Names must be unique")
	}
	return(.Call("checkHets", object))
}
.hetData <- setClass("hetData", contains = "list", validity = checkHets)