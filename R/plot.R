#' @include mpcross-class.R
#' @include geneticData-class.R
setMethod(f = "plot", signature = "mpcross", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"mpcross\". Use functions plotProbabilities or plotMosaic instead")
})
setMethod(f = "plot", signature = "geneticData", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"geneticData\". Use functions plotProbabilities or plotMosaic instead")
})
setMethod(f = "plot", signature = "probabilities", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"probabilities\". Use functions plotProbabilities or plotMosaic instead")
})
setMethod(f = "plot", signature = "imputed", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"imputed\". Use functions plotProbabilities or plotMosaic instead")
})

