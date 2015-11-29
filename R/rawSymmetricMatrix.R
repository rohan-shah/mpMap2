checkRawSymmetricMatrix <- function(object)
{
	errors <- c()
	if(any(object@levels > 0.5 | object@levels < 0))
	{
		errors <- c(errors, "Slot levels must contain values between 0 and 0.5")
	}
	if(length(object@data) != length(object@markers)*(length(object@markers)+1)/2)
	{
		browser()
		errors <- c(errors, "Slots markers and data had incompatible lengths")
	}
	if(any(diff(object@levels) < 0))
	{
		errors <- c(errors, "Values in slot levels must be in increasing order")
	}
	if(length(object@levels) >= 255)
	{
		errors <- c(errors, "At most 254 possible levels are allowed")
	}
	return(errors)
}
.rawSymmetricMatrix <- setClass("rawSymmetricMatrix", slots = list(data = "raw", markers = "character", levels = "numeric"), validity = checkRawSymmetricMatrix)
setMethod("[", signature(x = "rawSymmetricMatrix", i = "index", j = "index", drop = "logical"),
	function(x, i, j, ..., drop)
	{
		nMarkers <- length(x@markers)
		if(any((i > nMarkers) | (j > nMarkers) | (i < 1) | (j < 1))) stop("Indices were out of range")
		return(.Call("rawSymmetricMatrixSubsetIndices", x, i, j, drop, PACKAGE="mpMap2"))
	})
setMethod("[", signature(x = "rawSymmetricMatrix", i = "index", j = "index", drop = "missing"),
	function(x, i, j, ..., drop)
	{
		nMarkers <- length(x@markers)
		if(any((i > nMarkers) | (j > nMarkers) | (i < 1) | (j < 1))) stop("Indices were out of range")
		return(.Call("rawSymmetricMatrixSubsetIndices", x, i, j, TRUE, PACKAGE="mpMap2"))
	})
