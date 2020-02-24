#' @include mpcross-class.R
#' @include geneticData-class.R
#' @rdname redact
#' @title Redact sensitive information
#' 
#' This function redacts possibly sensitive information from objects, resulting in an object that is safe to publish. 
#'
#' @description Sensitive information includes names of genetic lines (both founding lines and final population lines) and marker names. All actual data (marker genotypes, imputed IBD genotypes, IBD probabilities, etc) are preserved. 
#' @param object The object of class \code{mpcross} to redact.
#' @return An object of class \code{mpcross}, with identifying information removed. 
#' @export
setGeneric(name = "redact", def = function(object){standardGeneric("redact")})
#' @rdname redact
setMethod(f = "redact", signature = "mpcross", definition = function(object)
{
	object@geneticData <- new("geneticDataList", lapply(object@geneticData, redact))
	return(object)
})
#' @rdname redact
setMethod(f = "redact", signature = "mpcrossRF", definition = function(object)
{
	redactedMpcross <- callNextMethod()
	return(new("mpcrossRF", redactedMpcross, rf = object@rf))
})
#' @rdname redact
setMethod(f = "redact", signature = "mpcrossLG", definition = function(object)
{
	redactedMpcross <- callNextMethod()
	return(new("mpcrossLG", redactedMpcross, lg = object@lg))
})
#' @rdname redact
setMethod(f = "redact", signature = "mpcrossMapped", definition = function(object)
{
	redactedMpcross <- callNextMethod()
	return(new("mpcrossMapped", redactedMpcross, rf = object@rf, map = object@map))
})
#' @rdname redact
setMethod(f = "redact", signature = "geneticData", definition = function(object)
{
	nFounders <- nFounders(object)
	originalFinalNames <- rownames(object@finals)
	originalFounderNames <- rownames(object@founders)
	nFinalLines <- nrow(object@finals)

	originalLineNames <- object@pedigree@lineNames
	newLineNames <- paste0("L", 1:length(object@pedigree@lineNames))
	#Keep the founder names
	newLineNames[1:nFounders] <- object@pedigree@lineNames[1:nFounders]
	names(newLineNames) <- originalLineNames

	rownames(object@finals) <- newLineNames[rownames(object@finals)]
	object@pedigree@lineNames <- as.vector(newLineNames)

	newFinalNames <- newLineNames[originalFinalNames]

	if(!is.null(object@imputed))
	{
		rownames(object@imputed@data) <- newLineNames[rownames(object@imputed@data)]
		if(!is.null(object@imputed@errors))
		{
			rownames(object@imputed@errors) <- newLineNames[rownames(object@imputed@errors)]
		}
	}
	if(!is.null(object@probabilities))
	{
		if(nrow(object@probabilities@key) == nFounders * nFounders)
		{
			nAlleles <- length(unique(object@probabilities@key[,3]))
			rownames(object@probabilities@data) <- paste0(rep(newFinalNames, each = nAlleles), " - ", rep(1:nAlleles, times = nFinalLines))
		}
		else
		{
			rownames(object@probabilities@data) <- paste0(rep(newFinalNames, each = nFounders), " - ", rep(originalFounderNames, times = nFinalLines))
		}
	}
	return(object)
})
