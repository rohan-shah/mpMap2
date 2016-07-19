#' @include pedigree-class.R
#' @include detailedPedigree-class.R
#' @export
setGeneric("selfing<-", function(object, value) standardGeneric("selfing<-"))
#' @export
setGeneric("selfing", function(object) standardGeneric("selfing"))
setMethod(f = "selfing", signature = "pedigree", definition = function(object)
{
	object@selfing
})
setReplaceMethod("selfing", "detailedPedigree", function(object, value)
{
	if(value != "infinite" && value != "finite")
	{
		stop("Selfing must be either \"finite\" or \"infinite\"")
	}
	object@selfing <- value
	object
})
setReplaceMethod("selfing", "pedigree", function(object, value)
{
	if(value != "infinite" && value != "finite")
	{
		stop("Selfing must be either \"finite\" or \"infinite\"")
	}
	object@selfing <- value
	object
})
