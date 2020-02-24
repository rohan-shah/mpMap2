checkRawSymmetricMatrix <- function(object)
{
	errors <- c()
	if(any(object@levels > 0.5 | object@levels < 0))
	{
		errors <- c(errors, "Slot levels must contain values between 0 and 0.5")
	}
	if(length(object@data) != length(object@markers)*(length(object@markers)+1)/2)
	{
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
	#Note that this creates a logical vector having the same length as object@data, before the any(...) is applied. Logicals are 4 bytes! So this is replaced with C code
	#if(any((object@data >= length(object@levels)) & object@data != as.raw(255)))
	if(.Call("checkRawSymmetricMatrix", object, PACKAGE="mpMap2"))
	{
		errors <- c(errors, "Value in slot data was too large")
	}
	return(errors)
}
.rawSymmetricMatrix <- setClass("rawSymmetricMatrix", slots = list(data = "raw", markers = "character", levels = "numeric"), validity = checkRawSymmetricMatrix)
#' @rdname internalOperators
#' @param i,j indices specifying elements to extract or replace
#' @param x object from which to extract element(s)
#' @param drop If TRUE the result is coerced to the lowest possible dimension
#' @param ... Currently unused
#' @return Various. Not for external use.
setMethod("[", signature(x = "rawSymmetricMatrix", i = "index", j = "index", drop = "logical"),
	function(x, i, j, ..., drop)
	{
		nMarkers <- length(x@markers)
		if(any(i > nMarkers) || any(j > nMarkers) || any(i < 1) || any(j < 1)) stop("Indices were out of range")
		return(.Call("rawSymmetricMatrixSubsetIndices", x, i, j, drop, PACKAGE="mpMap2"))
	})
#' @rdname internalOperators
setMethod("[", signature(x = "rawSymmetricMatrix", i = "index", j = "index", drop = "missing"),
	function(x, i, j, ..., drop)
	{
		nMarkers <- length(x@markers)
		if(any(i > nMarkers) || any(j > nMarkers) || any(i < 1) || any(j < 1)) stop("Indices were out of range")
		return(.Call("rawSymmetricMatrixSubsetIndices", x, i, j, TRUE, PACKAGE="mpMap2"))
	})
#' @rdname internalOperators
setMethod("[", signature(x = "rawSymmetricMatrix", i = "missing", j = "missing", drop = "missing"), 
	function(x, i, j, ..., drop)
	{
		data <- .Call("rawSymmetricMatrixUncompress", x, PACKAGE="mpMap2")
		rownames(data) <- colnames(data) <- x@markers
		return(data)
	})
#' @rdname internalOperators
setMethod("[", signature(x = "rawSymmetricMatrix", i = "matrix", j = "missing", drop = "missing"), 
	function(x, i, j, ..., drop)
	{
		if(ncol(i) != 2)
		{
			stop("Any matrix used for subsetting must have two columns")
		}
		nMarkers <- length(x@markers)
		if(any(i > nMarkers | i < 1))
		{
			stop("Indices were out of range")
		}
		if(!is.numeric(i))
		{
			stop("Any matrix used for subsetting must be numeric")
		}
		if(storage.mode(i) != "integer")
		{
			storage.mode(i) <- "integer"
		}
		return(.Call("rawSymmetricMatrixSubsetByMatrix", x, i, PACKAGE="mpMap2"))
	})
setAs("matrix", "rawSymmetricMatrix", def = function(from, to)
	{
		#This is mostly for testing purposes. Haven't bothered with writing efficient C code, because I don't think it can really be made that efficient. This function is mostly for testing purposes, giving us an easy way to generate objects of class rawSymmetricMatrix. 
		if(!isSymmetric(from))
		{
			stop("Input matrix must be symmetric")
		}
		if(is.null(rownames(from)) || is.null(colnames(from)) || !all.equal(rownames(from), colnames(from)))
		{
			stop("Row and column names are required, and must be identical")
		}
		levels <- sort(unique(as.vector(from)))
		values <- from[upper.tri(from, diag = TRUE)]
		indices <- match(values, levels)-1
		indices[is.na(values)] <- 0xff
		return(new("rawSymmetricMatrix", data = as.raw(indices), markers = rownames(from), levels = levels))
	})
setAs("rawSymmetricMatrix", "matrix", def = function(from, to)
	{
		return(from[])
	})
