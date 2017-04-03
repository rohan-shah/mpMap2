#' @include pedigree-class.R
#' @include detailedPedigree-class.R
#' @export
setGeneric("lineNames<-", function(object, value) standardGeneric("lineNames<-"))
#' @export
setGeneric("lineNames", function(object) standardGeneric("lineNames"))
setMethod(f = "lineNames", signature = "pedigree", definition = function(object)
{
	object@lineNames
})
setMethod(f = "lineNames", signature = "mpcross", definition = function(object)
{
	if(length(object@geneticData) == 1) return(rownames(finals(object)))
	else return(lapply(finals(object), rownames))
})
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

