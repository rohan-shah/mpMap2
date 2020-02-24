#' @include mpcross-class.R
#' @include geneticData-class.R
#' @title Plot methods
#' @description There are multiple meaningful ways to plot some mpMap2 objects. Please use \code{\link{plotProbabilities}} or \code{\link{plotMosaic}} instead. 
#' @rdname plot
#' @details There are multiple meaningful ways to plot some mpMap2 objects. In these cases the plot function is implemented but returns an error. Please use \code{\link{plotProbabilities}} or \code{\link{plotMosaic}} instead. 
#' @param x Unused
#' @param y Unused
#' @param ... Unused
#' @return None
setMethod(f = "plot", signature = "mpcross", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"mpcross\". Use functions plotProbabilities or plotMosaic instead")
})
#' @rdname plot
setMethod(f = "plot", signature = "geneticData", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"geneticData\". Use functions plotProbabilities or plotMosaic instead")
})
#' @rdname plot
setMethod(f = "plot", signature = "probabilities", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"probabilities\". Use functions plotProbabilities or plotMosaic instead")
})
#' @rdname plot
setMethod(f = "plot", signature = "imputed", definition = function(x, ...)
{
	stop("Function plot is not defined for an object of class \"imputed\". Use functions plotProbabilities or plotMosaic instead")
})

