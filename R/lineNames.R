#' @include pedigree-class.R
#' @include detailedPedigree-class.R
#' @title Get or set the genetic line names of a pedigree
#' @rdname lineNamesPedigree
#' @description Get or set the genetic line names of a pedigree
#' @details These functions get or set the names of the genetic lines associated with a pedigree. 
#' @export
setGeneric("lineNames<-", function(object, value) standardGeneric("lineNames<-"))
#' @rdname lineNamesGeneric
#' @title Get or set the genetic line names
#' @description Get or set the genetic line names associated with a pedigree or \code{mpcross} object.
#' @details These functions get or set the names of the genetic lines associated with a pedigree or \code{mpcross} object.
#' @export
setGeneric("lineNames", function(object) standardGeneric("lineNames"))
#' @rdname lineNamesPedigree
setMethod(f = "lineNames", signature = "pedigree", definition = function(object)
{
	object@lineNames
})
#' @include detailedPedigree-class.R
#' @title Get the genetic line names
#' @rdname lineNamesMpcross
#' @description Get the genetic line names of a population
#' @details These functions get the names of the genetic lines associated with an \code{mpcross} object.
#' @export
setMethod(f = "lineNames", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1) return(rownames(finals(object)))
	else return(lapply(finals(object), rownames))
})
#' @rdname lineNamesMpcross
setMethod(f = "lineNames", signature = "geneticData", definition = function(object)
{
	return(rownames(finals(object)))
})
#' @rdname lineNamesPedigree
setReplaceMethod("lineNames", "detailedPedigree", function(object, value)
{
	if(length(value) != length(object@lineNames))
	{
		stop("Input lineNames had the wrong length")
	}
	if(length(unique(value)) != length(value))
	{
		stop("Input lineNames cannot contain duplicates")
	}
	object@lineNames <- value
	object
})
#' @rdname lineNamesPedigree
setReplaceMethod("lineNames", "pedigree", function(object, value)
{
	if(length(value) != length(object@lineNames))
	{
		stop("Input lineNames had the wrong length")
	}
	if(length(unique(value)) != length(value))
	{
		stop("Input lineNames cannot contain duplicates")
	}
	object@lineNames <- value
	object
})

